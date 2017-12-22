/******************************************************************************
* Copyright (c) 2011, Michael P. Gerlek (mpg@flaxen.com)
*
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following
* conditions are met:
*
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in
*       the documentation and/or other materials provided
*       with the distribution.
*     * Neither the name of Hobu, Inc. or Flaxen Geo Consulting nor the
*       names of its contributors may be used to endorse or promote
*       products derived from this software without specific prior
*       written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
* COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
* BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
* OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
* AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
* OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
* OF SUCH DAMAGE.
****************************************************************************/

#pragma once

#include <list>

#include <pdal/Filter.hpp>
#include <pdal/Polygon.hpp>
#include <pdal/plugin.hpp>

#include "private/Point.hpp"

extern "C" int32_t CropFilter_ExitFunc();
extern "C" PF_ExitFunc CropFilter_InitPlugin();

namespace pdal
{

class ProgramArgs;
class GridPnp;

// removes any points outside of the given range
// updates the header accordingly
class PDAL_DLL CropFilter : public Filter
{
public:
    CropFilter();
    ~CropFilter();
    static void * create();
    static int32_t destroy(void *);
    std::string getName() const;

private:
    // This is just a way to marry a (multi)polygon with a list of its
    // GridPnp's that do the actual point-in-polygon operation.
    struct ViewGeom
    {
        ViewGeom(const Polygon& poly);
        ViewGeom(ViewGeom&& vg);

        Polygon m_poly;
        std::vector<std::unique_ptr<GridPnp>> m_gridPnps;
    };
    std::vector<Bounds> m_bounds;
    bool m_cropOutside;
    std::vector<Polygon> m_polys;
    SpatialReference m_assignedSrs;
    double m_distance;
    double m_distance2;
    std::vector<filter::Point> m_centers;
    std::vector<ViewGeom> m_geoms;
    std::vector<BOX2D> m_boxes;

    void addArgs(ProgramArgs& args);
    virtual void initialize();

    virtual void ready(PointTableRef table);
    virtual void spatialReferenceChanged(const SpatialReference& srs);
    virtual bool processOne(PointRef& point);
    virtual PointViewSet run(PointViewPtr view);
    bool crop(const PointRef& point, const BOX2D& box);
    void crop(const BOX2D& box, PointView& input, PointView& output);
    bool crop(const PointRef& point, GridPnp& g);
    void crop(const ViewGeom& g, PointView& input, PointView& output);
    bool crop(const PointRef& point, const filter::Point& center);
    void crop(const filter::Point& center, PointView& input,
        PointView& output);
    void transform(const SpatialReference& srs);

    CropFilter& operator=(const CropFilter&); // not implemented
    CropFilter(const CropFilter&); // not implemented
};

} // namespace pdal
