#pragma once

#include <random>
#include <vector>
#include <set>

#include "Comparison.hpp"
#include "Grid.hpp"
#include "VoxelRayTrace.hpp"

namespace pdal
{

struct grid_error : public std::runtime_error
{
    grid_error(const std::string& s) : std::runtime_error(s)
    {}
};

class GridPnp
{
public:
    using Point = std::pair<double, double>;

    // Initialize the point-in-poly engine by creating an overlay grid and
    // attaching the index of each polygon edge to any cell it passes
    // through.
    GridPnp(const std::vector<Point>& poly)
    {
        if (poly.size() < 4)
            throw grid_error("Invalid polygon: too few points.");
        if (poly[0] != poly[poly.size() - 1])
            throw grid_error("Polygon not closed.");

        m_poly.push_back(poly[0]);
        for (size_t i = 1; i < poly.size(); ++i)
        {
            const Point& p1 = poly[i];
            const Point& p2 = m_poly.back();

            if (Comparison::closeEnough(p1.first, p2.first) &&
                Comparison::closeEnough(p1.second, p2.second))
                continue;
            m_poly.push_back(p1);
        }
        if (m_poly.size() < 4)
            throw grid_error("Invalid polygon: too few points after "
                "culling consecutive duplicates.");

        // Also calculates the average edge length.
        double xAvgLen;
        double yAvgLen;
        calcBounds(xAvgLen, yAvgLen);
        XYIndex gridSize = calcGridSize(xAvgLen, yAvgLen);
        createGrid(gridSize);
        assignEdges();
    }


    bool inside(const Point& p)
    { return inside(p.first, p.second); }

    // Determine if a point is inside the polygon attached to this class.
    bool inside(double x, double y)
    {
        // Find the index of the grid cell at the position.
        // If the position isn't in the grid, we're certainly outside.
        XYIndex idx;
        if (!m_grid->cellPos(x, y, idx))
            return false;

        Cell& cell = m_grid->cell(idx);
        // If we don't already have a reference point with computed state,
        // do it now.
        if (!cell.computed())
            computeCell(cell, idx);

        // If there are no edges in the cell, the status of the cell is uniform,
        // just return the state.
        if (cell.empty())
            return cell.inside();
        return testCell(cell, x, y);
    }

private:
    using XYIndex = std::pair<size_t, size_t>;
    using Edge = std::pair<Point, Point>;

    enum class IntersectType
    {
        Cross,
        On,
        None
    };

    class Cell
    {
    public:
        Cell() : m_point(
            {std::numeric_limits<double>::quiet_NaN(),
             std::numeric_limits<double>::quiet_NaN() })
        {}

        void addEdge(size_t edge)
            { m_edges.push_back(edge); }
        bool empty() const
            { return m_edges.empty(); }
        void setPoint(double x, double y)
            { m_point = Point(x, y); }
        bool computed()
            { return !std::isnan(m_point.first); }
        GridPnp::Point point() const
            { return m_point; }
        const std::vector<size_t>& edges() const
            { return m_edges; }
        bool inside() const
            { return m_inside; }
        void setInside(bool inside)
            { m_inside = inside; }

    private:
        std::vector<size_t> m_edges;
        bool m_inside;
        GridPnp::Point m_point;
    };

    Point pointFromId(size_t id)
        { return Point{ m_poly[id].first, m_poly[id].second }; }
    Edge edgeFromId(size_t id)
        { return Edge {pointFromId(id), pointFromId(id + 1)}; }
    double xval(const Point& p)
        { return p.first; }
    double yval(const Point& p)
        { return p.second; }
    size_t numEdges()
        { return m_poly.size() - 1; }

    // Calculate the bounding box of the polygon.  At the same time
    // calculate the average length of the X and Y components of the
    // polygon edges.
    void calcBounds(double& xAvgLen, double& yAvgLen)
    {
        double xdist{0};
        double ydist{0};

        // Inialize max/min with X/Y of first point.
        Point p = pointFromId(0);
        m_xMin = xval(p);
        m_xMax = xval(p);
        m_yMin = yval(p);
        m_yMax = yval(p);

        // The first point is duplicated as the last, so we skip the last
        // point when looping.
        for (size_t id = 0; id < numEdges(); ++id)
        {
            Edge e = edgeFromId(id);
            Point& p1 = e.first;
            Point& p2 = e.second;

            // Calculate bounding box.
            m_xMin = std::min(m_xMin, xval(p1));
            m_xMax = std::max(m_xMax, xval(p1));
            m_yMin = std::min(m_yMin, yval(p1));
            m_yMax = std::max(m_yMax, yval(p1));

            // Sum the lengths of the X and Y components of the edges.
            xdist += std::abs(xval(p2) - xval(p1));
            ydist += std::abs(yval(p2) - yval(p1));
        }

        // Find the average X and Y component length.
        xAvgLen = xdist / numEdges();
        yAvgLen = ydist / numEdges();
    }


    // The paper calculates an X and Y based on a) the number of edges
    // and b) the relative length of edges in the X and Y direction.
    // This seems fine, but it misses out on considering the common
    // case where a polygon delineates some area of interest.  In this
    // case there is much "empty" space in the interior of the polygon and
    // it seems likely that most pnp tests will happen in the empty interior.
    // So, it would seem that we'd want the grid sufficiently large to cover
    // this case well.  That way, most pnp tests would take no work beyond
    // an array lookup since most cells would be empty.  Lots of tradeoff,
    // though, in preprocessing vs. actual pnp tests.  Hopefully more work
    // can be done on this later.  My stupid way of dealing with this is
    // to set a minimum grid size of 1000 cells.
    XYIndex calcGridSize(double xAvgLen, double yAvgLen)
    {
        // I'm setting a minmum number of cells as 1000, because, why not?
        size_t m = std::max(1000UL, m_poly.size());

        // See paper for this calc.
        double scalex = ((m_xMax - m_xMin) * yAvgLen) /
            ((m_yMax - m_yMin) * xAvgLen);
        double scaley = 1 / scalex;
        double mx = std::sqrt(m * scalex);
        double my = std::sqrt(m * scaley);

        // We always round up, because why not.
        return XYIndex(mx + 1, my + 1);
    }


    // Figure out the grid origin.
    void createGrid(XYIndex gridSize)
    {
        // Make the grid extend 1/2 cell beyond bounds box.
        double boxWidth = m_xMax - m_xMin;
        double boxHeight = m_yMax - m_yMin;
        //
        double cellWidth = boxWidth / (gridSize.first - 1);
        double cellHeight = boxHeight / (gridSize.second - 1);
        double xOrigin = m_xMin - (cellWidth / 2);
        double yOrigin = m_yMin - (cellHeight / 2);

        m_grid.reset(new Grid<Cell>(gridSize.first, gridSize.second,
            cellWidth, cellHeight, xOrigin, yOrigin));
        m_xDistribution.reset(
            new std::uniform_real_distribution<>(0, m_grid->cellWidth()));
        m_yDistribution.reset(
            new std::uniform_real_distribution<>(0, m_grid->cellHeight()));
    }


    // Loop through edges.  Add the edge to each cell traversed.
    void assignEdges()
    {
        for (size_t id = 0; id < numEdges(); ++id)
        {
            Edge e = edgeFromId(id);
            Point& p1 = e.first;
            Point& p2 = e.second;
            Point origin = m_grid->origin();
            VoxelRayTrace vrt(m_grid->cellWidth(), m_grid->cellHeight(),
                xval(origin), yval(origin),
                xval(p1), yval(p1), xval(p2), yval(p2));
            VoxelRayTrace::CellList traversedCells = vrt.emit();
            for (auto& c : traversedCells)
                m_grid->cell(XYIndex(c.first, c.second)).addEdge(id);
        }
    }


    // Determine if a point is collinear with an edge.
    bool pointCollinear(double x, double y, Edge edge)
    {
        Point p1 = edge.first;
        Point p2 = edge.second;
        double x1 = xval(p1);
        double x2 = xval(p2);
        double y1 = yval(p1);
        double y2 = yval(p2);

        // If p1 == p2, this will fail.

        // This is the same as saying slopes are equal.
        return Comparison::closeEnough((x - x2) * (y - y1),
            (y - y2) * (x - x1));
    }


    // Put a reference point in the cell.  Figure out if the reference point
    // is inside the polygon.
    void computeCell(Cell& cell, XYIndex& pos)
    {
        generateRefPoint(cell, pos);
        determinePointStatus(cell, pos);
    }


    // The paper uses point centers, but then has to deal with points that
    // are "singular" (rest on a polygon edge).  But there's nothing special
    // about the point center.  The center is just a point in a cell with
    // a known status (inside or outside the polygon).  So we just pick a point
    // that isn't collinear with any of the segments in the cell.  Eliminating
    // collinearity eliminates special cases when counting crossings.
    void generateRefPoint(Cell& cell, XYIndex& pos)
    {
        // A test point is valid if it's not collinear with any segments
        // in the cell.
        auto validTestPoint = [this](double x, double y, Cell& cell)
        {
            for (auto edge : cell.edges())
                if (pointCollinear(x, y, edgeFromId(edge)))
                    return false;
            return true;
        };

        Grid<Cell>::Point origin = m_grid->cellOrigin(pos);
        double x, y;
        do
        {
            x = xval(origin) + (*m_xDistribution)(m_ranGen);
            y = yval(origin) + (*m_yDistribution)(m_ranGen);
        } while (!validTestPoint(x, y, cell));
        cell.setPoint(x, y);
    }


    // Determine the status of a cell's reference point by drawing a segment
    // from the reference point in the cell to the left and count crossings.
    // Knowing the number of edge crossings and the inside/outside status
    // of the cell determines the status of this reference point.
    // If we're determining the status of the leftmost cell, choose a point
    // to the left of the leftmost cell, which is guaranteed to be outside
    // the polygon.
    void determinePointStatus(Cell& cell, XYIndex& pos)
    {
        Point p1(cell.point());

        double x1 = xval(p1);
        double y1 = yval(p1);

        size_t intersectCount = 0;
        if (pos.first == 0)
        {
            double x2 = x1 - m_grid->cellWidth();
            double y2 = y1;

            Edge edge{{x1, y1}, {x2, y2}};

            intersectCount = intersections(edge, cell.edges());
        }
        else
        {
            XYIndex prevPos {pos.first - 1, pos.second};
            Cell& prevCell = m_grid->cell(prevPos);
            if (!prevCell.computed())
                computeCell(prevCell, prevPos);
            double x2 = xval(prevCell.point());
            double y2 = yval(prevCell.point());

            Edge edge{{x1, y1}, {x2, y2}};

            // Stick the edges in the current cell and the previous cell
            // in a set so as not to double-count.
            std::set<size_t> edges;
            edges.insert(cell.edges().begin(), cell.edges().end());
            edges.insert(prevCell.edges().begin(), prevCell.edges().end());
            intersectCount = intersections(edge, edges);
            if (prevCell.inside())
                intersectCount++;
        }
        cell.setInside(intersectCount % 2 == 1);
    }


    // Determine the number of intersections between an edge and
    // all edges indexes by the 'edges' list.
    template<typename EDGES>
    size_t intersections(Edge& e1, const EDGES& edges)
    {
        size_t isect = 0;
        for (auto& edgeId : edges)
        {
            Edge e2 = edgeFromId(edgeId);
            if (intersects(e1, e2) != IntersectType::None)
                isect++;
        }
        return isect;
    }


    // Determine if a point in a cell is inside the polygon or outside.
    // We're always calling a point that lies on an edge as 'inside'
    // the polygon.
    bool testCell(Cell& cell, double x, double y)
    {
        Edge tester({x, y}, cell.point());

        bool inside = cell.inside();
        for (auto edgeIdx: cell.edges())
        {
            Edge other = edgeFromId(edgeIdx);
            IntersectType intersection = intersects(tester, other);
            if (intersection == IntersectType::On)
                return true;
            if (intersection == IntersectType::Cross)
                inside = !inside;
        }
        return inside;
    }

    // Determine if two edges interset.  Note that because of the way
    // we've chosen reference points, the two segments should never be
    // collinear, which eliminates some special cases.
    //
    // One segment endpoint lies on the other if the slope factor (t or u)
    // is one or 0 and the other factor is between 0 and 1.
    // This is standard math, but it's shown nicely on Stack Overflow
    // question 563198.  The variable names map to the good response there.
    IntersectType intersects(Edge& e1, Edge& e2)
    {
        using Vector = std::pair<double, double>;

        Vector p = e1.first;
        Vector r { e1.second.first - e1.first.first,
            e1.second.second - e1.first.second };
        Vector q = e2.first;
        Vector s { e2.second.first - e2.first.first,
            e2.second.second - e2.first.second };

        // Should never be 0.
        double rCrossS = r.first * s.second - r.second * s.first;
        Vector pq = { q.first - p.first, q.second - p.second };

        double pqCrossS = pq.first * s.second - pq.second * s.first;
        double t = (pqCrossS / rCrossS);
        bool tCloseEnough = Comparison::closeEnough(t, 0) ||
            Comparison::closeEnough(t, 1);
        bool intersect = (tCloseEnough || (t > 0 && t < 1));
        if (!intersect)
            return IntersectType::None;

        double pqCrossR = pq.first * r.second - pq.second * r.first;
        double u = (pqCrossR / rCrossS);
        bool uCloseEnough = Comparison::closeEnough(u, 0) ||
            Comparison::closeEnough(u, 1);
        intersect = (uCloseEnough || (u > 0 && u < 1));
        if (intersect)
        {
            if (uCloseEnough || tCloseEnough)
                return IntersectType::On;
            return IntersectType::Cross;
        }
        return IntersectType::None;
    }

    std::vector<Point> m_poly;
    std::mt19937 m_ranGen;
    std::unique_ptr<std::uniform_real_distribution<>> m_xDistribution;
    std::unique_ptr<std::uniform_real_distribution<>> m_yDistribution;
    std::unique_ptr<Grid<Cell>> m_grid;
    double m_xMin;
    double m_xMax;
    double m_yMin;
    double m_yMax;
};

} // namespace pdal
