
// Grid class that holds a
template <typename T>
class Grid
{
public:
    using Point = std::pair<double, double>;
    using Pos = std::pair<size_t, size_t>;

    // Could make all these ctor vars compile-time, but this probably
    // complicates things for users who don't have this information
    // at compiler time.
    Grid(size_t width, size_t height, double cellWidth, double cellHeight,
        double xOrigin = 0.0, double yOrigin = 0.0) :
            m_width(width), m_height(height),
            m_cellWidth(cellWidth), m_cellHeight(cellHeight),
            m_xOrigin(xOrigin), m_yOrigin(yOrigin), m_cells(width * height)
    {}
    T& cell(size_t xpos, size_t ypos)
        { return m_cells[index(xpos, ypos)]; }
    T& cell(Pos pos)
        { return cell(pos.first, pos.second); }

    // Convert an x/y location to a cell position.
    bool cellPos(double x, double y, Pos& pos)
    {
        x -= m_xOrigin;
        y -= m_yOrigin;
        if (x < 0 || y < 0)
            return false;
        pos.first = x / m_cellWidth;
        pos.second = y / m_cellHeight;
        return (pos.first < m_width && pos.second < m_height);
    }
    Point cellOrigin(size_t xpos, size_t ypos) const
    {
        return {m_xOrigin + m_cellWidth * xpos,
            m_yOrigin + m_cellHeight * ypos};
    }
    Point cellOrigin(Pos pos) const
        { return cellOrigin(pos.first, pos.second); }
    Point origin() const
        { return { m_xOrigin, m_yOrigin }; }
    double cellWidth() const
        { return m_cellWidth; }
    double cellHeight() const
        { return m_cellHeight; }

private:
    size_t m_width, m_height;
    double m_cellWidth, m_cellHeight;
    double m_xOrigin, m_yOrigin;
    std::vector<T> m_cells;

    size_t index(size_t xpos, size_t ypos)
        { return ypos * m_width + xpos; }
};
