/*

Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef PlanarDivisionRule_HPP_
#define PlanarDivisionRule_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCentreBasedDivisionRule.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "DomMeshBasedCellPopulationWithGhostNodes.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StromalCellMutationState.hpp"


// Forward declaration prevents circular include chain
// template<3, 3> class AbstractCentreBasedCellPopulation;
// template<3, 3> class AbstractCentreBasedDivisionRule;

/**
 * A class to generate two daughter cell positions, located a distance
 * AbstractCentreBasedCellPopulation::mMeinekeDivisionSeparation apart,
 * along a random axis. The midpoint between the two daughter cell
 * positions corresponds to the parent cell's position.
 */
// template<3, 3>
class PlanarDivisionRule : public AbstractCentreBasedDivisionRule<3, 3>
{
private:
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCentreBasedDivisionRule<3, 3> >(*this);
    }

    // DomMeshBasedCellPopulationWithGhostNodes<3>* mpStaticCastCellPopulation;

public:

    /**
     * Default constructor.
     */
    PlanarDivisionRule()
    {
    }

    /**
     * Empty destructor.
     */
    virtual ~PlanarDivisionRule()
    {
    }

    /**
     * Overridden CalculateCellDivisionVector() method.
     *
     * @param pParentCell  The cell to divide
     * @param rCellPopulation  The centre-based cell population
     *
     * @return the two daughter cell positions.
     */
    virtual std::pair<c_vector<double, 3>, c_vector<double, 3> > CalculateCellDivisionVector(CellPtr pParentCell,
        AbstractCentreBasedCellPopulation<3, 3>& rCellPopulation);
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(PlanarDivisionRule)

#endif // PlanarDivisionRule_HPP_
