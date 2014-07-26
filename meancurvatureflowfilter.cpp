/****************************************************************************
* MeshLab                                                           o o     *
* A versatile mesh processing toolbox                             o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *   
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#include "meancurvatureflowfilter.h"
#include <vcg/complex/algorithms/update/curvature_fitting.h>
#include <vcg/complex/algorithms/clean.h>
#include <QtScript>

MeanCurvaureFlowPlugin::MeanCurvaureFlowPlugin()
{ 
    typeList << FP_MOVE_VERTEX;

    foreach(FilterIDType tt , types())
        actionList << new QAction(filterName(tt), this);
}

QString MeanCurvaureFlowPlugin::filterName(FilterIDType filterId) const
{
    switch(filterId) {
    case FP_MOVE_VERTEX :  return QString("Mean curvature flow");
    default : assert(0);
    }
    return QString();
}

QString MeanCurvaureFlowPlugin::filterInfo(FilterIDType filterId) const
{
    switch(filterId) {
    case FP_MOVE_VERTEX :  return QString("Move the vertices of the mesh along vertex normal according to curvature.\n"
                                          " Warrning: This version works only with closed meshes.");
    default : assert(0);
    }
    return QString("Unknown Filter");
}

MeanCurvaureFlowPlugin::FilterClass MeanCurvaureFlowPlugin::getClass(QAction *a)
{
    switch(ID(a))
    {
    case FP_MOVE_VERTEX :  return MeshFilterInterface::Smoothing;
    default : assert(0);
    }
    return MeshFilterInterface::Generic;
}

void MeanCurvaureFlowPlugin::initParameterSet(QAction *action,MeshModel &m, RichParameterSet & parlst)
{
    switch(ID(action))
    {
        case FP_MOVE_VERTEX :
                                      parlst.addParam(new RichFloat("Time",
                                      0.01f,
                                      "Mean curvature factor",
                                      "Value which will be multiplied with curvature."));
        break;

        default : assert(0);
    }
    //Enable curvature and adjacency relations needed to compute a curvature.
    m.updateDataMask(MeshModel::MM_VERTCURV);
    m.updateDataMask(MeshModel::MM_VERTCURVDIR);
    m.updateDataMask(MeshModel::MM_VERTFACETOPO);
    m.updateDataMask(MeshModel::MM_FACEFACETOPO);

    vcg::tri::Clean<CMeshO>::RemoveUnreferencedVertex(m.cm);
    vcg::tri::Allocator<CMeshO>::CompactVertexVector(m.cm);
    vcg::tri::UpdateCurvatureFitting<CMeshO>::computeCurvature( m.cm );
}

// The Real Core Function doing the actual mesh processing.
bool MeanCurvaureFlowPlugin::applyFilter(QAction */*filter*/, MeshDocument &md, RichParameterSet & par, vcg::CallBackPos *cb)
{
    CMeshO &m = md.mm()->cm;

    for(unsigned int i = 0; i< m.vert.size(); i++)
    {
        m.vert[i].P() -= m.vert[i].cN() * (( m.vert[i].K1() + m.vert[i].K2() ) / 2.) * par.getFloat("Time");
    }
    // Log function dump textual info in the lower part of the MeshLab screen.
    Log("Successfully displaced %i vertices",m.vn);
    vcg::tri::UpdateBounding<CMeshO>::Box(m);
    return true;
}

MESHLAB_PLUGIN_NAME_EXPORTER(MeanCurvaureFlowPlugin)
