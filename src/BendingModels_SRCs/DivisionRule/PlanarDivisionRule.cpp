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

#include "PlanarDivisionRule.hpp"
#include "RandomNumberGenerator.hpp"
// #include "DomMeshBasedCellPopulationWithGhostNodes.hpp"
#include "Debug.hpp"

// template<3, 3>
std::pair<c_vector<double, 3>, c_vector<double, 3> > PlanarDivisionRule::CalculateCellDivisionVector(
    CellPtr pParentCell,
    AbstractCentreBasedCellPopulation<3, 3>& rCellPopulation)
{
    // mpStaticCastCellPopulation = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*>(&mrCellPopulation);


    // Get separation parameter
    // double separation = rCellPopulation.GetMeinekeDivisionSeparation();
    // double separation = SimulationTime::Instance()->GetTimeStep();
    double separation = 0.1;


    // Make a random direction vector of the required length
    c_vector<double, 3> random_vector;

    /*
     * Pick a random direction and move the parent cell backwards by 0.5*separation
     * in that direction and return the position of the daughter cell 0.5*separation
     * forwards in that direction.
     */
    // switch (SPACE_DIM)
    // {
        // case 1: //Random/Not implemented...
        // {
        //     double random_direction = -1.0 + 2.0*(RandomNumberGenerator::Instance()->ranf() < 0.5);

        //     random_vector(0) = 0.5*separation*random_direction;
        //     break;
        // }
        // case 2: //Random/Not implemented...
        // {
        //     double random_angle = 2.0*M_PI*RandomNumberGenerator::Instance()->ranf();

        //     random_vector(0) = 0.5*separation*cos(random_angle);
        //     random_vector(1) = 0.5*separation*sin(random_angle);
        //     break;
        // }
        // case 3:
        // {
            DomMeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*>(&rCellPopulation);

            unsigned parent_node_index = rCellPopulation.GetLocationIndexUsingCell(pParentCell);

            // Get pointer to this node
            Node<3>* p_node = rCellPopulation.GetNode(parent_node_index);

            // Loop over containing elements
            std::vector<unsigned> neighbouring_node_indices;
            
            for (typename Node<3>::ContainingElementIterator elem_iter = p_node->ContainingElementsBegin();
                elem_iter != p_node->ContainingElementsEnd();
                ++elem_iter)
            {
                // Get pointer to this containing element
                Element<3,3>* p_element = p_tissue->rGetMesh().GetElement(*elem_iter);
            
                // Loop over nodes contained in this element
                for (unsigned i=0; i<p_element->GetNumNodes(); i++)
                {
                    // Get index of this node and add its index to the set if not the original node
                    unsigned node_index = p_element->GetNodeGlobalIndex(i);

                    if(!p_tissue->IsGhostNode(node_index))
                    {
                        CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(node_index);

                        if (node_index != parent_node_index && p_cell->GetMutationState()->IsType<WildTypeCellMutationState>())
                        {
                            neighbouring_node_indices.push_back(node_index);
                        }
                    }
                }
            }
            // We do actually want the parent node in there too, after all!
            neighbouring_node_indices.push_back(parent_node_index);

            

            c_vector<double, 3> normal_vector;
	
            //	a00, a10, a20, a01, a11, a21, a02, a12, a22
            c_vector<double, 9> ATA;
            c_vector<double, 3> ATz;
            c_vector<double, 9> iATA;
            for (int ii=0; ii<9; ii++)
            {
                ATA[ii] = 0;
                iATA[ii] = 0;
                if(ii<3)
                {
                    ATz[ii] = 0;
                }
            }

            for (unsigned i=0; i<neighbouring_node_indices.size(); i++)
            {
                unsigned cell_i =neighbouring_node_indices[i];
                //int epithelialNodeIndex = second_order_neighs[i];
                // c_vector<double, 3> epithelial_location = mpExtendedMesh->GetNode(cell_i_ext)->rGetLocation();
                c_vector<double, 3> epithelial_location = (rCellPopulation.GetNode(cell_i))->rGetLocation();


                // Find Transpose(A)*A
                ATA[0] += epithelial_location[0]*epithelial_location[0];
                ATA[1] += epithelial_location[0]*epithelial_location[1];
                ATA[2] += epithelial_location[0];

                ATA[3] += epithelial_location[0]*epithelial_location[1];
                ATA[4] += epithelial_location[1]*epithelial_location[1];
                ATA[5] += epithelial_location[1];

                ATA[6] += epithelial_location[0];
                ATA[7] += epithelial_location[1];
                ATA[8] += 1.0;

                // Calculate Transpose(A)*z
                ATz[0] += epithelial_location[0]*epithelial_location[2];
                ATz[1] += epithelial_location[1]*epithelial_location[2];
                ATz[2] += epithelial_location[2];
            }
            
            // determinant of Transpose(A)*A
	        double ATA_det = ATA[0]*(ATA[4]*ATA[8]-ATA[5]*ATA[7]) - ATA[1]*(ATA[3]*ATA[8]-ATA[5]*ATA[6]) + ATA[2]*(ATA[3]*ATA[7]-ATA[6]*ATA[4]);

            // if ( (ATA_det >= pow(10,-10)) || (ATA_det <= -1.0*pow(10,-10)) )
            // {
                // Calculate the inverse of Transpose(A)*A
                iATA[0] = (1.0/ATA_det)*(ATA[4]*ATA[8] - ATA[7]*ATA[5]);
                iATA[1] = (1.0/ATA_det)*(ATA[2]*ATA[7] - ATA[8]*ATA[1]);
                iATA[2] = (1.0/ATA_det)*(ATA[1]*ATA[5] - ATA[4]*ATA[2]);

                iATA[3] = (1.0/ATA_det)*(ATA[5]*ATA[6] - ATA[8]*ATA[3]);
                iATA[4] = (1.0/ATA_det)*(ATA[0]*ATA[8] - ATA[6]*ATA[2]);
                iATA[5] = (1.0/ATA_det)*(ATA[2]*ATA[3] - ATA[5]*ATA[0]);

                iATA[6] = (1.0/ATA_det)*(ATA[3]*ATA[7] - ATA[6]*ATA[4]);
                iATA[7] = (1.0/ATA_det)*(ATA[1]*ATA[6] - ATA[7]*ATA[0]);
                iATA[8] = (1.0/ATA_det)*(ATA[0]*ATA[4] - ATA[3]*ATA[1]);
                
                // Calculate  normal = inverse(Transpose(A)*A)*(Transpose(A)*z)
                normal_vector[0] = iATA[0]*ATz[0] + iATA[1]*ATz[1] + iATA[2]*ATz[2];
                normal_vector[1] = iATA[3]*ATz[0] + iATA[4]*ATz[1] + iATA[5]*ATz[2];
                normal_vector[2] = iATA[6]*ATz[0] + iATA[7]*ATz[1] + iATA[8]*ATz[2];

                normal_vector = normal_vector/norm_2(normal_vector);
                //std::cout << "case 1" << "\n"; 

                // normal_vector[0] = -1.0*normal_vector[0];
                // normal_vector[1] = -1.0*normal_vector[1];
            // }

            random_vector[0] = 2.0*(RandomNumberGenerator::Instance()->ranf()) - 1.0;
            random_vector[1] = 2.0*(RandomNumberGenerator::Instance()->ranf()) - 1.0;
            random_vector[2] = 2.0*(RandomNumberGenerator::Instance()->ranf()) - 1.0;

            

            // PRINT_VECTOR(random_vector);
            // break;
        // }
    //     default:
    //         // This can't happen
    //         NEVER_REACHED;
    // }
    c_vector<double, 3> planar_vector;


    planar_vector = random_vector - (random_vector[0]*normal_vector[0] + random_vector[1]*normal_vector[1] + random_vector[2]*normal_vector[2])*normal_vector;
    
    planar_vector = separation*planar_vector/norm_2(planar_vector);
    
    // PRINT_VECTOR(planar_vector);
    
    c_vector<double, 3> parent_position = rCellPopulation.GetLocationOfCellCentre(pParentCell) - planar_vector;
    c_vector<double, 3> daughter_position = parent_position + planar_vector;

    // std::pair<c_vector<double, 3>, c_vector<double, 3> > positions(parent_position, daughter_position);
    std::pair<c_vector<double, 3>, c_vector<double, 3> > positions(daughter_position, parent_position);

    return positions;
}


// Explicit instantiation
// template class PlanarDivisionRule<1,1>;
// template class PlanarDivisionRule<1,2>;
// template class PlanarDivisionRule<2,2>;
// template class PlanarDivisionRule<1,3>;
// template class PlanarDivisionRule<2,3>;
// template class PlanarDivisionRule<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PlanarDivisionRule)
