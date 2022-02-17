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

#include "PeriodicRemeshCellsModifier.hpp"
#include "DomMeshBasedCellPopulationWithGhostNodes.hpp"
#include "VtkMeshWriter.hpp"

#include "AbstractCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StromalCellMutationState.hpp"
#include "Debug.hpp"



template<unsigned DIM>
PeriodicRemeshCellsModifier<DIM>::PeriodicRemeshCellsModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
    mOutputDirectory(""),
    mWidth(0.0),
    mDepth(0.0)
{
}

template<unsigned DIM>
void PeriodicRemeshCellsModifier<DIM>::SetOutputDirectory(std::string outputDirectory)
{
	mOutputDirectory = outputDirectory;
}

template<unsigned DIM>
void PeriodicRemeshCellsModifier<DIM>::SetWidth(double width)
{
	mWidth = width;
}

template<unsigned DIM>
void PeriodicRemeshCellsModifier<DIM>::SetDepth(double depth)
{
	mDepth = depth;
}

// std::string MeshModifier::GetOutputDirectory()
// {
// 	return mOutputDirectory;
// }

template<unsigned DIM>
PeriodicRemeshCellsModifier<DIM>::~PeriodicRemeshCellsModifier()
{
}

template<unsigned DIM>
void PeriodicRemeshCellsModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM>& rCellPopulation)
{
}

template<unsigned DIM>
void PeriodicRemeshCellsModifier<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM>& rCellPopulation)
{

    std::vector<Node<DIM>*> extended_nodes_no_ghosts(4*rCellPopulation.GetNumRealCells());
    // std::vector<Node<DIM>*> extended_nodes_no_ghosts(rCellPopulation.GetNumRealCells());
    int count = 0;

    MutableMesh<DIM,DIM>* mpMeshNoGhosts;
    if (mpMeshNoGhosts != NULL)
    {
    	delete mpMeshNoGhosts;
    }

    std::map<unsigned, unsigned> MeshNodeIndexMap;

    // Dom - Create a copy of original mesh
	for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // First, create and store a copy of this real node and cell

        unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        c_vector<double, DIM> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
		// PRINT_2_VARIABLES(real_node_index,count);

        // Create a copy of the node corresponding to this cell and store it
        Node<DIM>* p_real_node = new Node<DIM>(real_node_index, real_node_location);
        extended_nodes_no_ghosts[count] = p_real_node;
        MeshNodeIndexMap[count] = real_node_index;
        count++;
        
    }

    // Dom - Create a copy of original mesh
	for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // First, create and store a copy of this real node and cell

        unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        c_vector<double, DIM> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

        c_vector<double,DIM> image_node_location_width = real_node_location;
        if (real_node_location[0] >= mWidth*0.5)
        {
            image_node_location_width[0] -= mWidth;
        }
        else if (real_node_location[0] <  mWidth*0.5)
        {
            image_node_location_width[0] += mWidth;
        }

        Node<DIM>* p_real_node_width = new Node<DIM>(count, image_node_location_width);

        extended_nodes_no_ghosts[count] = p_real_node_width;
        MeshNodeIndexMap[count] = real_node_index;
        count++;


        c_vector<double,DIM> image_node_location_depth = real_node_location;
        if (real_node_location[1] >= mDepth*0.5)
        {
            image_node_location_depth[1] -= mDepth;
        }
        else if (real_node_location[1] <  mDepth*0.5)
        {
            image_node_location_depth[1] += mDepth;
        }

        Node<DIM>* p_real_node_depth = new Node<DIM>(count, image_node_location_depth);

        extended_nodes_no_ghosts[count] = p_real_node_depth;
        MeshNodeIndexMap[count] = real_node_index;
        count++;

        c_vector<double,DIM> image_node_location_edge = real_node_location;
        if (real_node_location[1] >= mDepth*0.5)
        {
            image_node_location_edge[1] -= mDepth;
        }
        else if (real_node_location[1] <  mDepth*0.5)
        {
            image_node_location_edge[1] += mDepth;
        }
        if (real_node_location[0] >= mWidth*0.5)
        {
            image_node_location_edge[0] -= mWidth;
        }
        else if (real_node_location[0] <  mWidth*0.5)
        {
            image_node_location_edge[0] += mWidth;
        }

        Node<DIM>* p_real_node_edge = new Node<DIM>(count, image_node_location_edge);

        extended_nodes_no_ghosts[count] = p_real_node_edge;
        MeshNodeIndexMap[count] = real_node_index;
        count++;
    }

    // for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
    //      cell_iter != rCellPopulation.End();
    //      ++cell_iter)
    // {
    //     unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
    //     c_vector<double, DIM> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

    //     c_vector<double,DIM> image_node_location_edge = real_node_location;
    //     if (real_node_location[1] >= mDepth*0.5)
    //     {
    //         image_node_location_edge[1] -= mDepth;
    //     }
    //     else if (real_node_location[1] <  mDepth*0.5)
    //     {
    //         image_node_location_edge[1] += mDepth;
    //     }
    //     if (real_node_location[0] >= mWidth*0.5)
    //     {
    //         image_node_location_edge[0] -= mWidth;
    //     }
    //     else if (real_node_location[0] <  mWidth*0.5)
    //     {
    //         image_node_location_edge[0] += mWidth;
    //     }

    //     Node<DIM>* p_real_node_edge = new Node<DIM>(count, image_node_location_edge);

    //     extended_nodes_no_ghosts[count] = p_real_node_edge;
    //     MeshNodeIndexMap[count] = real_node_index;
    //     count++;
    // }

    mpMeshNoGhosts = new MutableMesh<DIM,DIM>(extended_nodes_no_ghosts);

    // PRINT_2_VARIABLES(rCellPopulation.GetNumRealCells() , mpMeshNoGhosts->GetNumNodes());
    
    MutableMesh<2,3>* pMutableMesh = new MutableMesh<2,3>();

    std::vector<Node<3>*> nodes;

    DomMeshBasedCellPopulationWithGhostNodes<DIM>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);


    for (unsigned node_index=0; node_index<mpMeshNoGhosts->GetNumNodes(); node_index++)
    {
        const c_vector<double, 3>& r_location = mpMeshNoGhosts->GetNode(node_index)->rGetLocation();
        nodes.push_back(new Node<3>(node_index, r_location));
        pMutableMesh->AddNode(new Node<3>(node_index, r_location));
    }


    int element_iter = 0;
    double dist_threshold = 20.0;
    int num_nodes_original_mesh = rCellPopulation.GetNumRealCells();

    // PRINT_VARIABLE(mpMeshNoGhosts->GetNumElements());

    for (unsigned elem_index=0; elem_index<mpMeshNoGhosts->GetNumElements(); elem_index++) 
    { 
        bool element_contains_ghost_nodes = false;

        // Get a pointer to the element
        Element<DIM,DIM>* p_element = mpMeshNoGhosts->GetElement(elem_index);
        
        int Node0Index = p_element->GetNodeGlobalIndex(0);
        int Node1Index = p_element->GetNodeGlobalIndex(1);
        int Node2Index = p_element->GetNodeGlobalIndex(2);
        int Node3Index = p_element->GetNodeGlobalIndex(3);

        int node_intex[4] = {Node0Index, Node1Index, Node2Index, Node3Index};

        int node0GlobalIndex = MeshNodeIndexMap[Node0Index];
        int node1GlobalIndex = MeshNodeIndexMap[Node1Index];
        int node2GlobalIndex = MeshNodeIndexMap[Node2Index];
        int node3GlobalIndex = MeshNodeIndexMap[Node3Index];
        int node_global_intex[4] = {node0GlobalIndex, node1GlobalIndex, node2GlobalIndex, node3GlobalIndex};

        // PRINT_4_VARIABLES(Node0Index, Node1Index, Node2Index, Node3Index);
        // PRINT_2_VARIABLES(mpMeshNoGhosts->GetNumElements() , elem_index);

        // int number_of_ghosts = (p_tissue->IsGhostNode(node0GlobalIndex)) + (p_tissue->IsGhostNode(node1GlobalIndex)) + (p_tissue->IsGhostNode(node2GlobalIndex)) + (p_tissue->IsGhostNode(node3GlobalIndex));
        
        // int part_of_original_mesh = 0;
        // part_of_original_mesh = (Node0Index <= num_nodes_original_mesh) + (Node1Index <= num_nodes_original_mesh) + (Node2Index <= num_nodes_original_mesh) + (Node3Index <= num_nodes_original_mesh);
        // PRINT_VARIABLE(part_of_original_mesh);
        // if (part_of_original_mesh > 1 )
        if (true)
        {
            for(int j=0; j<4; j++)
            {
                CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(node_global_intex[(0+j)%4]);
                boost::shared_ptr<AbstractCellMutationState> p_state = p_cell->GetMutationState();

                bool is_0_stromal_cell = (p_state->IsType<StromalCellMutationState>()==true);

                if(is_0_stromal_cell) 
                // if(true) 
                {
                    CellPtr p_cell_1 = rCellPopulation.GetCellUsingLocationIndex(node_global_intex[(1+j)%4]);
                    boost::shared_ptr<AbstractCellMutationState> p_state_1 = p_cell_1->GetMutationState();

                    CellPtr p_cell_2 = rCellPopulation.GetCellUsingLocationIndex(node_global_intex[(2+j)%4]);
                    boost::shared_ptr<AbstractCellMutationState> p_state_2 = p_cell_2->GetMutationState();
                    
                    CellPtr p_cell_3 = rCellPopulation.GetCellUsingLocationIndex(node_global_intex[(3+j)%4]);
                    boost::shared_ptr<AbstractCellMutationState> p_state_3 = p_cell_3->GetMutationState();
                    
                    int number_of_epithelial_cell = (p_state_1->IsType<WildTypeCellMutationState>()==true) + (p_state_2->IsType<WildTypeCellMutationState>()==true) + (p_state_3->IsType<WildTypeCellMutationState>()==true);

                    if (number_of_epithelial_cell == 3)
                    // if (true)
                    {
                        // PRINT_3_VARIABLES(node_intex[(0+j)%4],node_intex[(1+j)%4],node_intex[(2+j)%4] );

                        // const c_vector<double, 3>& r_location1 = rCellPopulation.GetNode(node_global_intex[(1+j)%4])->rGetLocation();
                        // const c_vector<double, 3>& r_location2 = rCellPopulation.GetNode(node_global_intex[(2+j)%4])->rGetLocation();
                        // const c_vector<double, 3>& r_location3 = rCellPopulation.GetNode(node_global_intex[(3+j)%4])->rGetLocation();

                        // // c_vector<double, 3> difference12 = r_location1 - r_location2;
                        // double distance_between_nodes_12 = norm_2(r_location1 - r_location2);

                        // // c_vector<double, 3> difference13 = r_location1 - r_location3;
                        // double distance_between_nodes_13 = norm_2(r_location1 - r_location3);

                        // // c_vector<double, 3> difference23 = r_location2 - r_location3;
                        // double distance_between_nodes_23 = norm_2(r_location2 - r_location3);

                        // We do this because the nodes may move after, due to periodicity, but the tesselation happens at the begining of the timestep
                        // if(distance_between_nodes_12 <= dist_threshold && distance_between_nodes_13 <= dist_threshold && distance_between_nodes_23 <= dist_threshold)
                        if(true)
                        {
                            std::vector<Node<3>*> ElementNodes0;
                            ElementNodes0.push_back(nodes[node_intex[(1+j)%4]]);
                            ElementNodes0.push_back(nodes[node_intex[(2+j)%4]]);
                            ElementNodes0.push_back(nodes[node_intex[(3+j)%4]]);
                            pMutableMesh->AddElement(new Element<2,3>(element_iter,ElementNodes0));
                            element_iter++;

                        }
                    }
                    
                }
                
                
            }
        }


        
        

        // ITERATE OVER NODES owned by this element
        // for (unsigned local_index=0; local_index<4; local_index++)
        // {
        //     if (p_tissue->IsGhostNode(p_element->GetNodeGlobalIndex(local_index)) == true)
        //     {
        //         element_contains_ghost_nodes = true;
                
        //     }
        //     PRINT_VARIABLE(element_contains_ghost_nodes);
        // }
    }


    // for (unsigned node_index=0; node_index<rCellPopulation.rGetMesh().GetNumNodes(); node_index++)
    // {
    //     const c_vector<double, 3>& r_location = rCellPopulation.rGetMesh().GetNode(node_index)->rGetLocation();
    //     nodes.push_back(new Node<3>(node_index, r_location));
    //     pMutableMesh->AddNode(new Node<3>(node_index, r_location));
    // }

    // Test to make sure it does the thing we want - plot a single (random) triangle
    // std::vector<Node<3>*> ElementNodes;
    // ElementNodes.push_back(nodes[2]);
    // ElementNodes.push_back(nodes[5]);
    // ElementNodes.push_back(nodes[1]);
    // pMutableMesh->AddElement(new Element<2,3>(1,ElementNodes));


    // Can't copy the mesh over ... so lets try reconstruct it:
    // MutableMesh<DIM,DIM>* pOriginalMeshCopy  = new MutableMesh<DIM,DIM>();

    // MutableMesh<3,3>* pOriginalMesh = new MutableMesh<3,3>(nodes);
    // std::vector<Node<3>*> original_nodes;

    // for (unsigned node_index=0; node_index<pOriginalMesh->GetNumNodes(); node_index++)
    // {
    //     const c_vector<double, 3>& r_location = pOriginalMesh->GetNode(node_index)->rGetLocation();
    //     original_nodes.push_back(new Node<3>(node_index, r_location));
    // }

    // PRINT_2_VARIABLES(rCellPopulation.rGetMesh().GetNumNodes(), pOriginalMesh->GetNumNodes());

    // for (unsigned elem_index=0; elem_index<rCellPopulation.rGetMesh().GetNumElements(); elem_index++)
    // {
    //     Element<3,3>* pElement=rCellPopulation.rGetMesh().GetElement(elem_index);

    //     unsigned Node0Index = pElement->GetNodeGlobalIndex(0);
    //     unsigned Node1Index = pElement->GetNodeGlobalIndex(1);
    //     unsigned Node2Index = pElement->GetNodeGlobalIndex(2);
    //     unsigned Node3Index = pElement->GetNodeGlobalIndex(3);

    //     std::vector<Node<3>*> ElementNodes0;
    //     ElementNodes0.push_back(nodes[Node0Index]);
    //     ElementNodes0.push_back(nodes[Node1Index]);
    //     ElementNodes0.push_back(nodes[Node2Index]);
    //     ElementNodes0.push_back(nodes[Node3Index]);
    //     pOriginalMeshCopy->AddElement(new Element<3,3>(elem_index,ElementNodes0));

    // }
   		
        // int element_iter = 0;
        // double dist_threshold = 10.0;
        // for (unsigned elem_index=0; elem_index<pOriginalMesh->GetNumElements(); elem_index++)
        // {   
        //     // PRINT_VARIABLE(elem_index);

        //     Element<3,3>* pElement=pOriginalMesh->GetElement(elem_index);

        //     // PRINT_VARIABLE(pElement);

        //     for( int j=0; j<1; j++) //for( int j=0; j<4; j++)
        //     {
        //         // unsigned Node0Index = pElement->GetNodeGlobalIndex((0+j)%4);
        //         // unsigned Node1Index = pElement->GetNodeGlobalIndex((1+j)%4);
        //         // unsigned Node2Index = pElement->GetNodeGlobalIndex((2+j)%4);
        //         // unsigned Node3Index = pElement->GetNodeGlobalIndex((3+j)%4);

        //         unsigned Node0Index = pElement->GetNodeGlobalIndex(0);
        //         unsigned Node1Index = pElement->GetNodeGlobalIndex(1);
        //         unsigned Node2Index = pElement->GetNodeGlobalIndex(2);
        //         unsigned Node3Index = pElement->GetNodeGlobalIndex(3);

        //         // PRINT_4_VARIABLES(Node0Index,Node1Index,Node2Index,Node3Index);

        //         int number_of_ghosts = 0;//(rCellPopulation.rGetMesh().IsGhostNode(Node0Index)) + (rCellPopulation.rGetMesh().IsGhostNode(Node1Index)) + (rCellPopulation.rGetMesh().IsGhostNode(Node2Index)) + (rCellPopulation.rGetMesh().IsGhostNode(Node3Index));
        //         // PRINT_VARIABLE(number_of_ghosts);

        //         if (number_of_ghosts == 0)
        //         {

        //             const c_vector<double, 3>& r_location0 = pOriginalMesh->GetNode(Node0Index)->rGetLocation();
        //             const c_vector<double, 3>& r_location1 = pOriginalMesh->GetNode(Node1Index)->rGetLocation();
        //             const c_vector<double, 3>& r_location2 = pOriginalMesh->GetNode(Node2Index)->rGetLocation();
        //             const c_vector<double, 3>& r_location3 = pOriginalMesh->GetNode(Node3Index)->rGetLocation();


        //             c_vector<double, 3> difference01 = pOriginalMesh->GetVectorFromAtoB(r_location0, r_location1);
        //             double distance_between_nodes_01 = norm_2(difference01);

        //             c_vector<double, 3> difference02 = pOriginalMesh->GetVectorFromAtoB(r_location0, r_location2);
        //             double distance_between_nodes_02 = norm_2(difference02);

        //             c_vector<double, 3> difference03 = pOriginalMesh->GetVectorFromAtoB(r_location0, r_location3);
        //             double distance_between_nodes_03 = norm_2(difference03);

        //             c_vector<double, 3> difference12 = pOriginalMesh->GetVectorFromAtoB(r_location1, r_location2);
        //             double distance_between_nodes_12 = norm_2(difference12);

        //             c_vector<double, 3> difference13 = pOriginalMesh->GetVectorFromAtoB(r_location1, r_location3);
        //             double distance_between_nodes_13 = norm_2(difference13);

        //             c_vector<double, 3> difference23 = pOriginalMesh->GetVectorFromAtoB(r_location2, r_location3);
        //             double distance_between_nodes_23 = norm_2(difference23);


        //             if(distance_between_nodes_01 <= dist_threshold && distance_between_nodes_03 <= dist_threshold && distance_between_nodes_13 <= dist_threshold)
        //             {
        //                 std::vector<Node<3>*> ElementNodes0;
        //                 ElementNodes0.push_back(original_nodes[Node0Index]);
        //                 ElementNodes0.push_back(original_nodes[Node1Index]);
        //                 ElementNodes0.push_back(original_nodes[Node3Index]);
        //                 pMutableMesh->AddElement(new Element<2,3>(element_iter,ElementNodes0));
        //                 element_iter++;
        //             }

        //             if(distance_between_nodes_01 <= dist_threshold && distance_between_nodes_02 <= dist_threshold && distance_between_nodes_12 <= dist_threshold)
        //             {
        //                 std::vector<Node<3>*> ElementNodes1;
        //                 ElementNodes1.push_back(original_nodes[Node0Index]);
        //                 ElementNodes1.push_back(original_nodes[Node1Index]);
        //                 ElementNodes1.push_back(original_nodes[Node2Index]);
        //                 pMutableMesh->AddElement(new Element<2,3>(element_iter,ElementNodes1));
        //                 element_iter++;
        //             }

        //             if(distance_between_nodes_02 <= dist_threshold && distance_between_nodes_03 <= dist_threshold && distance_between_nodes_23 <= dist_threshold)
        //             {
        //                 std::vector<Node<3>*> ElementNodes2;
        //                 ElementNodes2.push_back(original_nodes[Node0Index]);
        //                 ElementNodes2.push_back(original_nodes[Node2Index]);
        //                 ElementNodes2.push_back(original_nodes[Node3Index]);
        //                 pMutableMesh->AddElement(new Element<2,3>(element_iter,ElementNodes2));
        //                 element_iter++;
        //             }

        //             if(distance_between_nodes_12 <= dist_threshold && distance_between_nodes_13 <= dist_threshold && distance_between_nodes_23 <= dist_threshold)
        //             {
        //                 std::vector<Node<3>*> ElementNodes3;
        //                 ElementNodes3.push_back(original_nodes[Node1Index]);
        //                 ElementNodes3.push_back(original_nodes[Node2Index]);
        //                 ElementNodes3.push_back(original_nodes[Node3Index]);
        //                 pMutableMesh->AddElement(new Element<2,3>(element_iter,ElementNodes3));
        //                 element_iter++;
        //             }
        //         }
                

        //         // CellPtr p_cell = this->GetCellUsingLocationIndex(node_index);

                

        //         // PRINT_4_VARIABLES(Node0Index, Node1Index, Node2Index, Node3Index);

        //         // if ( all nodes are epithlial )
        //         // {   
                
        //         //     std::vector(Node<3>*) ElementNodes;
            
        //         //     unsigned node_elem_indices(3) = {Node0Index Node1Index Node2Index};
                
        //         //     for i=1:3
        //         //     {
        //         //         ElementNodes.push_back(pMutableMesh->GetNode(node_elem_indices(i));
        //         //     }
        //         //     pMutableMesh->AddElement(new Element<2,3>(element_iter, ElementNodes));
        //         //     element_iter++;
        //         // }

        //     }
        // }



    // int element_iter = 0; 
    // for (unsigned elem_index=0; elem_index<rCellPopulation.rGetMesh().GetNumElements(); elem_index++)
    // {    

    //     Element<3,3>* pElement=rCellPopulation.rGetMesh().GetElement(elem_index);

    //     PRINT_VARIABLE(pElement);

    //     // for j=0:3
    //     // {
    //     //     Node0Index = pElement->GetNodeGLobalIndex((0+j)%3);
    //     //     Node1Index = pElement->GetNodeGLobalIndex((1+j)%3);
    //     //     Node2Index = pElement->GetNodeGLobalIndex((2+j)%3);

    //     //     // if ( all nodes are epithlial )
    //     //     // {   
            
    //     //     //     std::vector(Node<3>*) ElementNodes;
        
    //     //     //     unsigned node_elem_indices(3) = {Node0Index Node1Index Node2Index};
            
    //     //     //     for i=1:3
    //     //     //     {
    //     //     //         ElementNodes.push_back(pMutableMesh->GetNode(node_elem_indices(i));
    //     //     //     }
    //     //     //     pMutableMesh->AddElement(new Element<2,3>(element_iter, ElementNodes));
    //     //     //     element_iter++;
    //     //     // }
    //     // }
    // }


    // int element_iter = 0; 
    // for (unsigned node_index=0; node_index<rCellPopulation.rGetMesh().GetNumNodes(); node_index++)
    // {
    //     // PRINT_VARIABLE(rCellPopulation.GetNode(node_index));

    //     Node<DIM>* p_node = rCellPopulation.GetNode(node_index);

    //     // Get the cell corresponding to this node
    //     CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(node_index);

    //     // Get mutation state of cell
    // 	boost::shared_ptr<AbstractCellMutationState> p_state = p_cell->GetMutationState();

    // 	// Get whether cell is dead
    // 	bool cell_is_dead = p_cell->IsDead();

    // 	// Get whether this cell is a live epithelial cell
    // 	bool is_live_epithelial_cell = (p_state->IsType<WildTypeCellMutationState>()==true) && !cell_is_dead;
		
	// 	DomMeshBasedCellPopulationWithGhostNodes<DIM>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

    // 	if (is_live_epithelial_cell)
    // 	{
    // 	    // Iterate over elements of Mesh containing this node and no ghost nodes
	// 		c_vector<unsigned, 3> tri_el;

	// 		// This needs both stromal and epithelial cells... so cant do monolayer..
    // 		for (typename Node<DIM>::ContainingElementIterator elem_iter = p_node->ContainingElementsBegin();
	// 	         elem_iter != p_node->ContainingElementsEnd();
	// 	         ++elem_iter)
    // 		{
	// 			// Get a pointer to the element
	// 			Element<DIM,DIM>* p_element = rCellPopulation.rGetMesh().GetElement(*elem_iter);
				
	// 			// ITERATE OVER NODES owned by this element
	// 			std::vector<unsigned> temp_triangular_element;
	// 			bool is_element_connected_to_tissue = false;
	// 			for (unsigned local_index=0; local_index<4; local_index++)
	// 			{
	// 				unsigned neighbour_index = p_element->GetNodeGlobalIndex(local_index);

	// 				CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(neighbour_index);

	// 				if (p_cell->GetMutationState()->IsType<WildTypeCellMutationState>())
	// 				{
	// 					temp_triangular_element.push_back(neighbour_index);
	// 				}
	// 				else if(p_cell->GetMutationState()->IsType<StromalCellMutationState>() || p_tissue->IsGhostNode(neighbour_index))
	// 				{
	// 					is_element_connected_to_tissue = true;
	// 				}

	// 			}

	// 			// This means, we only consider adding if theres 3 epithelial cells AND 1 stromal cell
	// 			if(temp_triangular_element.size() == 3 && is_element_connected_to_tissue)
	// 			{
	// 				std::sort(temp_triangular_element.begin(), temp_triangular_element.end());
	// 				// tri_el[0] = temp_triangular_element[0];
	// 				// tri_el[1] = temp_triangular_element[1];
	// 				// tri_el[2] = temp_triangular_element[2];

    //                 std::vector<Node<3>*> ElementNodes;
    //                 ElementNodes.push_back(nodes[temp_triangular_element[0]]);
    //                 ElementNodes.push_back(nodes[temp_triangular_element[1]]);
    //                 ElementNodes.push_back(nodes[temp_triangular_element[2]]);
    //                 pMutableMesh->AddElement(new Element<2,3>(element_iter, ElementNodes));
    //                 element_iter++;



	// 				// if(put_element_in_tri == 1)
	// 				// {
	// 				// 	epithelial_triangulation.push_back(tri_el);
	// 				// }

	// 			}
	// 		}
			
	// 	}
    // }











    std::ostringstream time_string;
    time_string << SimulationTime::Instance()->GetTimeStepsElapsed();
    std::string results_file = "periodic_mesh_results_" + time_string.str();
    VtkMeshWriter<2,3>* p_vtk_mesh_writer = new VtkMeshWriter<2,3>(mOutputDirectory, results_file, false);
    p_vtk_mesh_writer->WriteFilesUsingMesh(*pMutableMesh);

    // delete pOriginalMesh;
}



template<unsigned DIM>
void PeriodicRemeshCellsModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateAtEndOfOutputTimeStep(rCellPopulation);


}


template<unsigned DIM>
void PeriodicRemeshCellsModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
// template class MeshModifier<1>;
// template class MeshModifier<2>;
template class PeriodicRemeshCellsModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PeriodicRemeshCellsModifier)

