#include "DensityDependantCellKillerStrip3DWithGhostNodesSimpleV2.hpp"
#include "Debug.hpp"

/* Apoptosis for cells that are epithelial that go outside the box by the cutoff length from the centre:
 *	 ___________________________________
 *	|	   	  |				  |		    |
 *	|		  |				  |		    |
 *	|  Death  |<---(x0,y0)--->|  Death  |
 *	|		  |			 	  |		    |
 *	|_________|_______________|_________|
 */
DensityDependantCellKillerStrip3DWithGhostNodesSimpleV2::DensityDependantCellKillerStrip3DWithGhostNodesSimpleV2(AbstractCellPopulation<3>* pCrypt, double cut_off, double density_threshold, double cellPopulationWidth, double cellPopulationDepth)
    : AbstractCellKiller<3>(pCrypt),
    mCellsRemovedByAnoikis(0),
	mCutOffLength(cut_off),
	mDensityThreshold(density_threshold),
	mCellPopulationWidth(cellPopulationWidth),
    mCellPopulationDepth(cellPopulationDepth)
{
    // Sets up output file
//	OutputFileHandler output_file_handler(mOutputDirectory + "AnoikisData/", false);
//	mAnoikisOutputFile = output_file_handler.OpenOutputFile("results.anoikis");
}

DensityDependantCellKillerStrip3DWithGhostNodesSimpleV2::~DensityDependantCellKillerStrip3DWithGhostNodesSimpleV2()
{
//    mAnoikisOutputFile->close();
}

void DensityDependantCellKillerStrip3DWithGhostNodesSimpleV2::SetOutputDirectory(std::string outputDirectory)
{
	mOutputDirectory = outputDirectory;

	OutputFileHandler output_file_handler(mOutputDirectory+"/", false);
    out_stream deathLocationFile = output_file_handler.OpenOutputFile("densityDeaths.dat");
    *deathLocationFile << "time \t";
    for (unsigned i=0; i<3; i++)
    {
        *deathLocationFile << "location" << i << "\t";
    }
    *deathLocationFile << "Cell Density " << "\t";
    *deathLocationFile << "Cell ID " << "\n";
    deathLocationFile->close();
}

std::string DensityDependantCellKillerStrip3DWithGhostNodesSimpleV2::GetOutputDirectory()
{
	return mOutputDirectory;
}

/*
 * Method to get the neighbouring nodes (excluding ghost nodes) of a particular node
 * Can then be used to identify the type of cells that surround a particular cell.
 */
std::set<unsigned> DensityDependantCellKillerStrip3DWithGhostNodesSimpleV2::GetNeighbouringNodeIndices(MutableMesh<3,3>* ExtendedMesh, std::map<unsigned, unsigned> ExtendedMeshNodeIndexMap, unsigned nodeIndex)
{
	// Create a set of neighbouring node indices
	std::set<unsigned> neighbouring_node_indices;

	// Need access to the mesh but can't get to it because the cell killer only owns a
	// pointer to an AbstractCellPopulation
    DomMeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*> (this->mpCellPopulation);

	// Find the indices of the elements owned by this node
    std::set<unsigned> containing_elem_indices = ExtendedMesh->GetNode(nodeIndex)->rGetContainingElementIndices();

	// CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(nodeIndex);
	// double x_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[0];
	// double y_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[1];
	// double z_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[2];

	//double cut_off = 1.5;

    // Iterate over these elements
    for (std::set<unsigned>::iterator elem_iter=containing_elem_indices.begin();
         elem_iter != containing_elem_indices.end();
	     ++elem_iter)
    {
	    // Get all the nodes contained in this element
		unsigned neighbour_extended_index;
	    unsigned neighbour_global_index;

	    for (unsigned local_index=0; local_index<4; local_index++)
	    {
	    	neighbour_extended_index = ExtendedMesh->GetElement(*elem_iter)->GetNodeGlobalIndex(local_index);
			neighbour_global_index = ExtendedMeshNodeIndexMap[neighbour_extended_index];

	    	// Don't want to include the original node or ghost nodes
			if(!p_tissue->IsGhostNode(neighbour_global_index))
			{
				// CellPtr p_cell_n = p_tissue->GetCellUsingLocationIndex(neighbour_global_index);
				// double x_location_n = this->mpCellPopulation->GetLocationOfCellCentre(p_cell_n)[0];
				// double y_location_n = this->mpCellPopulation->GetLocationOfCellCentre(p_cell_n)[1];
				// double z_location_n = this->mpCellPopulation->GetLocationOfCellCentre(p_cell_n)[2];

				// bool less_than_cutt_off = pow((x_location-x_location_n),2)+pow((y_location-y_location_n),2)+pow((z_location-z_location_n),2) < pow(mCutOffLength,2);

				if( (neighbour_global_index != nodeIndex))
				{
					neighbouring_node_indices.insert(neighbour_extended_index);
				}
			}
			// else if(p_tissue->IsGhostNode(neighbour_global_index))
			// {
			// 	//const c_vector<double, 3>& node_location = this->GetNode(neighbour_global_index)->rGetLocation();
			// 	//bool less_than_cutt_off = sqrt(pow((x_location-node_location[0]),2)+pow((y_location-node_location[1]),2)+pow((z_location-node_location[2]),2)) < mCutOffLength;

			// 	if( (neighbour_global_index != nodeIndex))// && less_than_cutt_off)
			// 	{
			// 		neighbouring_node_indices.insert(neighbour_global_index);
			// 	}
			// }
	    }
    }
    return neighbouring_node_indices;		// This will contain repeats, but it doesn't matter
}

/** Method to determine if an epithelial cell has lost all contacts with the stromal cells below
 * TRUE if cell has popped up
 * FALSE if cell remains in the monolayer
 */
c_vector<double,2> DensityDependantCellKillerStrip3DWithGhostNodesSimpleV2::IsCellTooSmall(MutableMesh<3,3>* ExtendedMesh, std::map<unsigned, unsigned> ExtendedMeshNodeIndexMap, unsigned nodeIndex, std::vector<unsigned> first_order_neighs)
{
	c_vector<double,2> is_too_small_density;

	DomMeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*> (this->mpCellPopulation);

	bool has_cell_popped_up = false;	// Initialising
	is_too_small_density[0] = 0.0;

   	unsigned num_stromal_neighbours = 0;
	unsigned num_ghost_neighbours = 0;
	unsigned num_epithelial_neighbours = 0;
	double average_cell_distance = 0.0;
	
	bool has_neighbours_begun_apoptosis = false;

	CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(nodeIndex);
	double cell_area_p = p_cell->GetCellData()->GetItem("cell_area");

	double x_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[0];
	double y_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[1];
	double z_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[2];

   	// Iterate over the neighbouring cells to check the number of stromal cell neighbours
	for(unsigned j=0; j<first_order_neighs.size(); j++)
	{
		unsigned neighbour_iter = first_order_neighs[j];
		unsigned global_index = ExtendedMeshNodeIndexMap[neighbour_iter];
		CellPtr p_cell_n = p_tissue->GetCellUsingLocationIndex(global_index);

		// We dont want cells to die if its neighbours are dying
		if(p_cell_n->HasApoptosisBegun())
		{
			has_neighbours_begun_apoptosis = true;
		}

		// Check how many stromal cells its conected to
		if ( (!p_tissue->IsGhostNode(global_index))
				&& (p_tissue->GetCellUsingLocationIndex(global_index)->GetMutationState()->IsType<StromalCellMutationState>()==true) )
   		{
			num_stromal_neighbours += 1;
		}
		// else if ( (!p_tissue->IsGhostNode(global_index))
		// 		&& (p_tissue->GetCellUsingLocationIndex(global_index)->GetMutationState()->IsType<WildTypeCellMutationState>()==true)
		// 		&& !(p_cell_n->HasApoptosisBegun()) )
		else if ( (!p_tissue->IsGhostNode(global_index))
				&& (p_tissue->GetCellUsingLocationIndex(global_index)->GetMutationState()->IsType<WildTypeCellMutationState>()==true))
   		{
     		c_vector<double, 3> epithelial_location = ExtendedMesh->GetNode(neighbour_iter)->rGetLocation();

			average_cell_distance = average_cell_distance + sqrt( pow(x_location - epithelial_location[0],2) + pow(y_location - epithelial_location[1],2) + pow(z_location - epithelial_location[2],2));
			num_epithelial_neighbours += 1;
		}
		// Check how many ghost nodes its conected to
		else if (p_tissue->IsGhostNode(global_index))
		{
			num_ghost_neighbours += 1;
		}
   	}

	average_cell_distance = (average_cell_distance/num_epithelial_neighbours);
	is_too_small_density[1] = cell_area_p;

	// Add cell density in cell data here.

   	if( (cell_area_p < mDensityThreshold) && (has_neighbours_begun_apoptosis == false))
   	{
   		has_cell_popped_up = true;
		   is_too_small_density[0] = 1.0;
   	}
	//PRINT_VARIABLE(num_stromal_neighbours);
   
    // PRINT_2_VARIABLES(average_cell_distance,has_cell_popped_up);
	// return has_cell_popped_up;
	return is_too_small_density;
}

std::vector<c_vector<unsigned, 3> > DensityDependantCellKillerStrip3DWithGhostNodesSimpleV2::GetEpithelialMesh(MutableMesh<3,3>* ExtendedMesh, std::map<unsigned, unsigned> ExtendedMeshNodeIndexMap)
{
    // Get a pointer to this node in mpExtendedMesh
    std::vector<c_vector<unsigned, 3> > epithelial_triangulation;

	DomMeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*> (this->mpCellPopulation);
	for (unsigned elem_index=0; elem_index<ExtendedMesh->GetNumElements(); elem_index++) 
	{ 

		// Get a pointer to the element
		Element<3,3>* p_element = ExtendedMesh->GetElement(elem_index);
			
		unsigned Node0Index = p_element->GetNodeGlobalIndex(0);
		unsigned Node1Index = p_element->GetNodeGlobalIndex(1);
		unsigned Node2Index = p_element->GetNodeGlobalIndex(2);
		unsigned Node3Index = p_element->GetNodeGlobalIndex(3);
		unsigned node_index[4] = {Node0Index, Node1Index, Node2Index, Node3Index};

		unsigned node0GlobalIndex = ExtendedMeshNodeIndexMap[Node0Index];
		unsigned node1GlobalIndex = ExtendedMeshNodeIndexMap[Node1Index];
		unsigned node2GlobalIndex = ExtendedMeshNodeIndexMap[Node2Index];
		unsigned node3GlobalIndex = ExtendedMeshNodeIndexMap[Node3Index];
		unsigned node_global_intex[4] = {node0GlobalIndex, node1GlobalIndex, node2GlobalIndex, node3GlobalIndex}; //check for 4 unique

		unsigned number_of_ghosts = (p_tissue->IsGhostNode(node0GlobalIndex)) + (p_tissue->IsGhostNode(node1GlobalIndex)) + (p_tissue->IsGhostNode(node2GlobalIndex)) + (p_tissue->IsGhostNode(node3GlobalIndex));

		if(number_of_ghosts == 0)
		{
			for(unsigned j=0; j<4; j++)
			{
				c_vector<unsigned, 3> tri_el = zero_vector<double>(3);

				CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(node_global_intex[(0+j)%4]);
				boost::shared_ptr<AbstractCellMutationState> p_state = p_cell->GetMutationState();

				bool is_0_stromal_cell = (p_state->IsType<StromalCellMutationState>()==true);

				if(is_0_stromal_cell)
				{
					CellPtr p_cell_1 = p_tissue->GetCellUsingLocationIndex(node_global_intex[(1+j)%4]);
					boost::shared_ptr<AbstractCellMutationState> p_state_1 = p_cell_1->GetMutationState();

					CellPtr p_cell_2 = p_tissue->GetCellUsingLocationIndex(node_global_intex[(2+j)%4]);
					boost::shared_ptr<AbstractCellMutationState> p_state_2 = p_cell_2->GetMutationState();
							
					CellPtr p_cell_3 = p_tissue->GetCellUsingLocationIndex(node_global_intex[(3+j)%4]);
					boost::shared_ptr<AbstractCellMutationState> p_state_3 = p_cell_3->GetMutationState();
							
					unsigned number_of_epithelial_cell = (p_state_1->IsType<WildTypeCellMutationState>()==true) + (p_state_2->IsType<WildTypeCellMutationState>()==true) + (p_state_3->IsType<WildTypeCellMutationState>()==true);

					if (number_of_epithelial_cell == 3)
					{
						tri_el[0] = node_index[(1+j)%4];
						tri_el[1] = node_index[(2+j)%4];
						tri_el[2] = node_index[(3+j)%4];

						epithelial_triangulation.push_back(tri_el);
					}
							
				}
						
			}
		}

	}

	return epithelial_triangulation;
}

std::vector<c_vector<unsigned, 10> > DensityDependantCellKillerStrip3DWithGhostNodesSimpleV2::GetEpithelialNeighbours(std::vector<c_vector<unsigned, 3> > rEpithelialMeshVector, unsigned number_of_cells)
{
	std::vector<c_vector<unsigned, 10> > node_neighbours;

	// c_vector<unsigned, 10> tmp_neighbours = zero_vector<unsigned>(10);

	unsigned temp_neighbours[4*number_of_cells][10];
	for(unsigned i=0; i<4*number_of_cells; i++)
	{
		for(unsigned j=0; j<10; j++)
		{
			temp_neighbours[i][j] = 0;
		}
	}

	for(unsigned i=0; i<rEpithelialMeshVector.size(); i++)
	{
		unsigned node_A = rEpithelialMeshVector[i][0];
		unsigned node_B = rEpithelialMeshVector[i][1];
		unsigned node_C = rEpithelialMeshVector[i][2];

		// Do node_A first
		unsigned row_neighs_A[10];
		for(unsigned j=0; j<10; j++)
		{
			row_neighs_A[j] = temp_neighbours[node_A][j];
		}
		if(std::find(std::begin(row_neighs_A), std::end(row_neighs_A), node_B) == std::end(row_neighs_A))
		{
			for (unsigned iter_A = 0; iter_A<10; iter_A++) 
			{
				if(temp_neighbours[node_A][iter_A] == 0)
				{
					temp_neighbours[node_A][iter_A] = node_B;
					break;
				}
			}
		}
		for(unsigned j=0; j<10; j++)
		{
			row_neighs_A[j] = temp_neighbours[node_A][j];
		}
		if(std::find(std::begin(row_neighs_A), std::end(row_neighs_A), node_C) == std::end(row_neighs_A))
		{
			for (unsigned iter_A = 0; iter_A<10; iter_A++) 
			{
				if(temp_neighbours[node_A][iter_A] == 0)
				{
					temp_neighbours[node_A][iter_A] = node_C;
					break;
				}
			}
		}


		// Then do node_B
		unsigned row_neighs_B[10];
		for(unsigned j=0; j<10; j++)
		{
			row_neighs_B[j] = temp_neighbours[node_B][j];
		}
		if(std::find(std::begin(row_neighs_B), std::end(row_neighs_B), node_A) == std::end(row_neighs_B))
		{
			for (unsigned iter_B = 0; iter_B<10; iter_B++) 
			{
				if(temp_neighbours[node_B][iter_B] == 0)
				{
					temp_neighbours[node_B][iter_B] = node_A;
					break;
				}
			}
		}
		for(unsigned j=0; j<10; j++)
		{
			row_neighs_B[j] = temp_neighbours[node_B][j];
		}
		if(std::find(std::begin(row_neighs_B), std::end(row_neighs_B), node_C) == std::end(row_neighs_B))
		{
			for (unsigned iter_B = 0; iter_B<10; iter_B++) 
			{
				if(temp_neighbours[node_B][iter_B] == 0)
				{
					temp_neighbours[node_B][iter_B] = node_C;
					break;
				}
			}
		}

		// Then do node_C
		unsigned row_neighs_C[10];
		for(unsigned j=0; j<10; j++)
		{
			row_neighs_C[j] = temp_neighbours[node_C][j];
		}
		if(std::find(std::begin(row_neighs_C), std::end(row_neighs_C), node_A) == std::end(row_neighs_C))
		{
			for (unsigned iter_C = 0; iter_C<10; iter_C++) 
			{
				if(temp_neighbours[node_C][iter_C] == 0)
				{
					temp_neighbours[node_C][iter_C] = node_A;
					break;
				}
			}
		}
		for(unsigned j=0; j<10; j++)
		{
			row_neighs_C[j] = temp_neighbours[node_C][j];
		}
		if(std::find(std::begin(row_neighs_C), std::end(row_neighs_C), node_B) == std::end(row_neighs_C))
		{
			for (unsigned iter_C = 0; iter_C<10; iter_C++) 
			{
				if(temp_neighbours[node_C][iter_C] == 0)
				{
					temp_neighbours[node_C][iter_C] = node_B;
					break;
				}
			}
		}

	}
	
	for(unsigned i=0; i<4*number_of_cells; i++)
	{
		std::vector<unsigned> holder_neighbour;
		c_vector<unsigned, 10> another_holder_neighbour = zero_vector<double>(10);
		for(unsigned j=0; j<10; j++)
		{
			holder_neighbour.push_back(temp_neighbours[i][j]);

		}
		
		another_holder_neighbour[0] = holder_neighbour[0];
		another_holder_neighbour[1] = holder_neighbour[1];
		another_holder_neighbour[2] = holder_neighbour[2];
		another_holder_neighbour[3] = holder_neighbour[3];
		another_holder_neighbour[4] = holder_neighbour[4];
		another_holder_neighbour[5] = holder_neighbour[5];
		another_holder_neighbour[6] = holder_neighbour[6];
		another_holder_neighbour[7] = holder_neighbour[7];
		another_holder_neighbour[8] = holder_neighbour[8];
		another_holder_neighbour[9] = holder_neighbour[9];

		node_neighbours.push_back(another_holder_neighbour);

	}
	return node_neighbours;

}

/** A method to return a vector that indicates which cells should be killed by anoikis
 */
std::vector<c_vector<unsigned,2> > DensityDependantCellKillerStrip3DWithGhostNodesSimpleV2::RemoveByAnoikis()
{
	std::vector<c_vector<unsigned,2> > cells_to_remove;

	return cells_to_remove;
}

/* Cell Killer that kills transit cells that move beyond the walls of the box
 * and also any transit cells that pop upwards and become detached from the stromal
 * cells
*/
void DensityDependantCellKillerStrip3DWithGhostNodesSimpleV2::CheckAndLabelCellsForApoptosisOrDeath()
{
    // Get the information at this timestep for each node index that says whether to remove by anoikis 
    DomMeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*> (this->mpCellPopulation);
	double domain_tollerance = 0.1;


	// unsigned num_cells = p_tissue->GetNumRealCells();
    // std::vector<Node<3>*> extended_nodes(4*num_cells);
	// std::map<unsigned, unsigned> mExtendedMeshNodeIndexMap;
	// MutableMesh<3,3>* mpExtendedMesh = nullptr;

	// unsigned count = 0;
	// // Dom - Create a copy of original mesh
    // for (typename AbstractCellPopulation<3>::Iterator cell_iter = mpCellPopulation->Begin();
    //     cell_iter != mpCellPopulation->End();
    //     ++cell_iter)
    // {
    //     // First, create and store a copy of this real node and cell
    //     unsigned real_node_index = mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
    //     c_vector<double, 3> real_node_location = mpCellPopulation->GetLocationOfCellCentre(*cell_iter);

    //     // Create a copy of the node corresponding to this cell and store it
    //     Node<3>* p_real_node = new Node<3>(real_node_index, real_node_location);
    //     extended_nodes[count] = p_real_node;

    //     // Populate mExtendedMeshNodeIndexMap
    //     mExtendedMeshNodeIndexMap[count] = real_node_index;

    //     count++;
    // }

	unsigned num_cells = p_tissue->GetNumRealCells();
	unsigned num_nodes = mpCellPopulation->GetNumNodes();
    std::vector<Node<3>*> extended_nodes(4*num_nodes);
	std::map<unsigned, unsigned> mExtendedMeshNodeIndexMap;
	MutableMesh<3,3>* mpExtendedMesh = nullptr;

	unsigned count = 0;
	// Dom - Create a copy of original mesh
    for (unsigned i=0; i<num_nodes; i++)
    {
        // First, create and store a copy of this real node and cell
        unsigned real_node_index = i;
        c_vector<double, 3> real_node_location = mpCellPopulation->GetNode(real_node_index)->rGetLocation();

        // Create a copy of the node corresponding to this cell and store it
        Node<3>* p_real_node = new Node<3>(real_node_index, real_node_location);
        extended_nodes[count] = p_real_node;

        // Populate mExtendedMeshNodeIndexMap
        mExtendedMeshNodeIndexMap[count] = real_node_index;

        count++;
    }


    // First, extend the mesh in the x-direction
    // for (typename AbstractCellPopulation<3>::Iterator cell_iter = mpCellPopulation->Begin();
    //     cell_iter != mpCellPopulation->End();
    //     ++cell_iter)
    // {
    //     // First, create and store a copy of this real node and cell
    //     unsigned real_node_index = mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
    //     c_vector<double, 3> real_node_location = mpCellPopulation->GetLocationOfCellCentre(*cell_iter);

    //     // Create a copy of the node corresponding to this cell and store it
    //      Node<3>* p_real_node = new Node<3>(real_node_index, real_node_location);

    //     // Compute the location of the image node corresponding to this node
    //     c_vector<double,3> image_node_location = real_node_location;
    //     if (real_node_location[0] >= mCellPopulationWidth*0.5)
    //     {
    //         image_node_location[0] -= mCellPopulationWidth;
    //     }
    //     else if (real_node_location[0] <  mCellPopulationWidth*0.5)
    //     {
    //         image_node_location[0] += mCellPopulationWidth;
    //     }

    //     // Create a copy of the node corresponding to this cell, suitable translated, and store it
    //     Node<3>* p_image_node = new Node<3>(count, image_node_location);
    //     extended_nodes[count] = p_image_node;

    //     // Populate mExtendedMeshNodeIndexMap
    //     mExtendedMeshNodeIndexMap[count] = real_node_index;

    //     count++;
    // }
	for (unsigned i=0; i<num_nodes; i++)
    {
        // First, create and store a copy of this real node and cell
        unsigned real_node_index = i;
        c_vector<double, 3> real_node_location = mpCellPopulation->GetNode(real_node_index)->rGetLocation();

        // Compute the location of the image node corresponding to this node
        c_vector<double,3> image_node_location = real_node_location;
        if (real_node_location[0] >= mCellPopulationWidth*0.5)
        {
            image_node_location[0] -= mCellPopulationWidth;
        }
        else if (real_node_location[0] <  mCellPopulationWidth*0.5)
        {
            image_node_location[0] += mCellPopulationWidth;
        }

        // Create a copy of the node corresponding to this cell, suitable translated, and store it
        Node<3>* p_image_node = new Node<3>(count, image_node_location);
        extended_nodes[count] = p_image_node;

        // Populate mExtendedMeshNodeIndexMap
        mExtendedMeshNodeIndexMap[count] = real_node_index;

        count++;
    }


    // Second, extend this extended mesh in the y-direction
    // (We don't need to store the real nodes anymore)
    // for (typename AbstractCellPopulation<3>::Iterator cell_iter = mpCellPopulation->Begin();
    //     cell_iter != mpCellPopulation->End();
    //     ++cell_iter)
    // {
    //     // First, create and store a copy of this real node and cell
    //     unsigned real_node_index = mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
    //     c_vector<double, 3> real_node_location = mpCellPopulation->GetLocationOfCellCentre(*cell_iter);

    //     // Compute the location of the image node corresponding to this node
    //     c_vector<double,3> image_node_location = real_node_location;

    //     if (real_node_location[1] >= mCellPopulationDepth*0.5)
    //     {
    //         image_node_location[1] -= mCellPopulationDepth;
    //     }
    //     else if (real_node_location[1] <  mCellPopulationDepth*0.5)
    //     {
    //         image_node_location[1] += mCellPopulationDepth;
    //     }

    //     // Create a copy of the node corresponding to this cell, suitable translated, and store it
    //     Node<3>* p_image_node = new Node<3>(count, image_node_location);
    //     extended_nodes[count] = p_image_node;

    //     // Populate mExtendedMeshNodeIndexMap
    //     mExtendedMeshNodeIndexMap[count] = real_node_index;

    //     count++;
    // }
	for (unsigned i=0; i<num_nodes; i++)
    {
        // First, create and store a copy of this real node and cell
        unsigned real_node_index = i;
        c_vector<double, 3> real_node_location = mpCellPopulation->GetNode(real_node_index)->rGetLocation();

        // Compute the location of the image node corresponding to this node
        c_vector<double,3> image_node_location = real_node_location;

        if (real_node_location[1] >= mCellPopulationDepth*0.5)
        {
            image_node_location[1] -= mCellPopulationDepth;
        }
        else if (real_node_location[1] <  mCellPopulationDepth*0.5)
        {
            image_node_location[1] += mCellPopulationDepth;
        }

        // Create a copy of the node corresponding to this cell, suitable translated, and store it
        Node<3>* p_image_node = new Node<3>(count, image_node_location);
        extended_nodes[count] = p_image_node;

        // Populate mExtendedMeshNodeIndexMap
        mExtendedMeshNodeIndexMap[count] = real_node_index;

        count++;
    }

    // Thirdly, extend this extended mesh so that we cover the corners too
    // (We don't need to store the real nodes anymore)
    // for (typename AbstractCellPopulation<3>::Iterator cell_iter = mpCellPopulation->Begin();
    //     cell_iter != mpCellPopulation->End();
    //     ++cell_iter)
    // {
    //     // First, create and store a copy of this real node and cell
    //     unsigned real_node_index = mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
    //     c_vector<double, 3> real_node_location = mpCellPopulation->GetLocationOfCellCentre(*cell_iter);

    //     // Compute the location of the image node corresponding to this node
    //     c_vector<double,3> image_node_location = real_node_location;

    //     if (real_node_location[1] >= mCellPopulationDepth*0.5)
    //     {
    //         image_node_location[1] -= mCellPopulationDepth;
    //     }
    //     else if (real_node_location[1] <  mCellPopulationDepth*0.5)
    //     {
    //         image_node_location[1] += mCellPopulationDepth;
    //     }
	// 	if (real_node_location[0] >= mCellPopulationWidth*0.5)
    //     {
    //         image_node_location[0] -= mCellPopulationWidth;
    //     }
    //     else if (real_node_location[0] <  mCellPopulationWidth*0.5)
    //     {
    //         image_node_location[0] += mCellPopulationWidth;
    //     }

    //     // Create a copy of the node corresponding to this cell, suitable translated, and store it
    //     Node<3>* p_image_node = new Node<3>(count, image_node_location);
    //     extended_nodes[count] = p_image_node;

    //     // Populate mExtendedMeshNodeIndexMap
    //     mExtendedMeshNodeIndexMap[count] = real_node_index;

    //     count++;
    // }
	for (unsigned i=0; i<num_nodes; i++)
    {
        // First, create and store a copy of this real node and cell
        unsigned real_node_index = i;
        c_vector<double, 3> real_node_location = mpCellPopulation->GetNode(real_node_index)->rGetLocation();

        // Compute the location of the image node corresponding to this node
        c_vector<double,3> image_node_location = real_node_location;

        if (real_node_location[1] >= mCellPopulationDepth*0.5)
        {
            image_node_location[1] -= mCellPopulationDepth;
        }
        else if (real_node_location[1] <  mCellPopulationDepth*0.5)
        {
            image_node_location[1] += mCellPopulationDepth;
        }
		if (real_node_location[0] >= mCellPopulationWidth*0.5)
        {
            image_node_location[0] -= mCellPopulationWidth;
        }
        else if (real_node_location[0] <  mCellPopulationWidth*0.5)
        {
            image_node_location[0] += mCellPopulationWidth;
        }

        // Create a copy of the node corresponding to this cell, suitable translated, and store it
        Node<3>* p_image_node = new Node<3>(count, image_node_location);
        extended_nodes[count] = p_image_node;

        // Populate mExtendedMeshNodeIndexMap
        mExtendedMeshNodeIndexMap[count] = real_node_index;

        count++;
    }

    // We now construct mpExtendedMesh using extended_nodes
    if (mpExtendedMesh != NULL)
    {
    	delete mpExtendedMesh;
    }
    mpExtendedMesh = new MutableMesh<3,3>(extended_nodes);
	// TRACE("Extended Mesh");
	//This needs fixing
	//   assert(p_tissue->GetVoronoiTessellation()!=NULL);	// This fails during archiving of a simulation as Voronoi stuff not archived yet

    std::vector<c_vector<unsigned,2> > cells_to_remove;
    c_vector<unsigned,2> individual_node_information;	// Will store the node index and whether to remove or not (1 or 0)

	std::vector<c_vector<unsigned, 3> > epithelial_triangulation = GetEpithelialMesh(mpExtendedMesh, mExtendedMeshNodeIndexMap);
	// TRACE("Got Ep Mesh");
	// std::vector<c_vector<unsigned, 10> > epithelial_neighbours = GetEpithelialNeighbours(epithelial_triangulation, num_cells);
	std::vector<c_vector<unsigned, 10> > epithelial_neighbours = GetEpithelialNeighbours(epithelial_triangulation, num_nodes);
	// TRACE("Got Neighbours");

	// PRINT_VARIABLE(SimulationTime::Instance()->GetTime());

	// for (AbstractCellPopulation<3>::Iterator cell_iter = p_tissue->Begin();
    // 	 cell_iter != p_tissue->End();
    // 	 ++cell_iter)
	for (typename AbstractCellPopulation<3>::Iterator cell_iter = mpCellPopulation->Begin();
        cell_iter != mpCellPopulation->End();
        ++cell_iter)
	{

		unsigned node_index = mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);


		// unsigned node_index = p_tissue->GetNodeCorrespondingToCell(*cell_iter)->GetIndex();
		//assert((!p_tissue->IsGhostNode(node_index)));
		if(!p_tissue->IsGhostNode(node_index))
		{
			// Initialise
			individual_node_information[0] = node_index;
			individual_node_information[1] = 0;


			// CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(node_index);
			CellPtr p_cell = mpCellPopulation->GetCellUsingLocationIndex(node_index);


			// if (cell_iter->GetMutationState()->IsType<StromalCellMutationState>()==false && cell_iter->GetMutationState()->IsType<WildTypeCellMutationState>()==true)
			if (p_cell->GetMutationState()->IsType<StromalCellMutationState>()==false && p_cell->GetMutationState()->IsType<WildTypeCellMutationState>()==true)
			{
				// assert(cell_iter->GetMutationState()->IsType<WildTypeCellMutationState>()==true);

				double x_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[0];
				double y_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[1];
				double z_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[2];


				////////////////////////////////////////////////////////////////////////////////////////////////

					// Find the first order neighbours

					// std::vector<unsigned> first_order_neighs_vect;
					// for(unsigned j=0; j<10; j++)
					// {
					// 	first_order_neighs_vect.push_back(epithelial_neighbours[node_index][j]);
					// 	// first_order_neighs_vect.push_back(epithelial_neighbours[cell_i_ext][j]);
					// }
					// // order and remove doubles
					// std::sort(first_order_neighs_vect.begin(), first_order_neighs_vect.end());
					// first_order_neighs_vect.erase(std::unique(first_order_neighs_vect.begin(), first_order_neighs_vect.end()), first_order_neighs_vect.end());
					// // Remove zero from vector
					// first_order_neighs_vect.erase(std::remove(first_order_neighs_vect.begin(), first_order_neighs_vect.end(), 0), first_order_neighs_vect.end());
					// first_order_neighs_vect.erase(std::remove(first_order_neighs_vect.begin(), first_order_neighs_vect.end(), node_index), first_order_neighs_vect.end());

					// // Determining whether to remove this cell of it's too small

					// c_vector<double,2> is_too_small_density = this->IsCellTooSmall(mpExtendedMesh,mExtendedMeshNodeIndexMap,node_index, first_order_neighs_vect);

					// PRINT_3_VARIABLES(node_index, is_too_small_density[1], first_order_neighs_vect.size() );
					// PRINT_VECTOR(first_order_neighs_vect);

				////////////////////////////////////////////////////////////////////////////////////////////////


				// Do killing in the x-direction
				// if ( pow(x_location - 0.5*mCellPopulationWidth,2) >= pow(mCutOffLength,2) && !(cell_iter->IsDead()) && (p_cell->GetAge() > 1) )
				if ( pow(x_location - 0.5*mCellPopulationWidth,2) >= pow(mCutOffLength,2) && !(p_cell->IsDead()) && (p_cell->GetAge() >= 1) )
				{
					// // Find the first order neighbours

					std::vector<unsigned> first_order_neighs_vect;
					for(unsigned j=0; j<10; j++)
					{
						first_order_neighs_vect.push_back(epithelial_neighbours[node_index][j]);
						// first_order_neighs_vect.push_back(epithelial_neighbours[cell_i_ext][j]);
					}
					// order and remove doubles
					std::sort(first_order_neighs_vect.begin(), first_order_neighs_vect.end());
					first_order_neighs_vect.erase(std::unique(first_order_neighs_vect.begin(), first_order_neighs_vect.end()), first_order_neighs_vect.end());
					// Remove zero from vector
					first_order_neighs_vect.erase(std::remove(first_order_neighs_vect.begin(), first_order_neighs_vect.end(), 0), first_order_neighs_vect.end());
					first_order_neighs_vect.erase(std::remove(first_order_neighs_vect.begin(), first_order_neighs_vect.end(), node_index), first_order_neighs_vect.end());

					// Determining whether to remove this cell of it's too small

					c_vector<double,2> is_too_small_density = this->IsCellTooSmall(mpExtendedMesh,mExtendedMeshNodeIndexMap,node_index, first_order_neighs_vect);
					
					// PRINT_3_VARIABLES(node_index, is_too_small_density[1], first_order_neighs_vect.size() );
					
					// if( (this->IsCellTooSmall(mpExtendedMesh,mExtendedMeshNodeIndexMap,node_index, first_order_neighs_vect)==true) )
										
					if( is_too_small_density[0] == 1.0)
					{
						if( !(p_cell->HasApoptosisBegun()) )
						{
							p_cell->StartApoptosis();
							individual_node_information[1] = 1;

							SimulationTime* p_time = SimulationTime::Instance();
							c_vector<double, 3> cell_location = p_tissue->GetNode(node_index)->rGetLocation();

							OutputFileHandler output_file_handler(mOutputDirectory+"/", false);
							out_stream deathLocationFile = output_file_handler.OpenOutputFile("densityDeaths.dat", std::ios::app);

							*deathLocationFile << p_time->GetTime() << "\t";
							for (unsigned i=0; i<3; i++)
							{
								*deathLocationFile << cell_location[i] << "\t";
							}
							*deathLocationFile << is_too_small_density[1] << "\t";
							*deathLocationFile << node_index << "\n";
							deathLocationFile->close();

						}

					}

				}

			}

			cells_to_remove.push_back(individual_node_information);
		}
	}

	delete mpExtendedMesh;

    // Keep a record of how many cells have been removed at this timestep
    this->SetNumberCellsRemovedByAnoikis(cells_to_remove);
    this->SetLocationsOfCellsRemovedByAnoikis(cells_to_remove);

    // Need to avoid trying to kill any cells twice (i.e. both by anoikis or random apoptosis)
    // Loop over these vectors individually and kill any cells that they tell you to

    // for (unsigned i=0; i<cells_to_remove.size(); i++)
    // {
    // 	if (cells_to_remove[i][1] == 1)
    // 	{
	// 		// PRINT_VARIABLE("cell removed");
    // 		// Get cell associated to this node
			
    // 		CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(cells_to_remove[i][0]);
	// 		if( !(p_cell->HasApoptosisBegun()) )
	// 		{
	// 			p_cell->StartApoptosis();

	// 		}

    // 	}
    // }

	for (AbstractCellPopulation<3>::Iterator cell_iter = p_tissue->Begin();
    	 cell_iter != p_tissue->End();
    	 ++cell_iter)
	{
		unsigned node_index = p_tissue->GetNodeCorrespondingToCell(*cell_iter)->GetIndex();
		CellPtr p_cell_A = mpCellPopulation->GetCellUsingLocationIndex(node_index);

		if(!p_tissue->IsGhostNode(node_index) && !(p_cell_A->IsDead()))
		{

			if(p_cell_A->HasApoptosisBegun())
			{
				// PRINT_2_VARIABLES(cell_iter->GetTimeUntilDeath(),SimulationTime::Instance()->GetTime());
				if (p_cell_A->GetTimeUntilDeath() <= 2*SimulationTime::Instance()->GetTimeStep())
				{
					p_cell_A->Kill();
					// TRACE("Cell Removed By Density");

					// SimulationTime* p_time = SimulationTime::Instance();
            		// c_vector<double, 3> cell_location = p_tissue->GetNode(node_index)->rGetLocation();

					// OutputFileHandler output_file_handler(mOutputDirectory+"/", false);
					// out_stream deathLocationFile = output_file_handler.OpenOutputFile("densityDeaths.dat", std::ios::app);

					// *deathLocationFile << p_time->GetTime() << "\t";
					// for (unsigned i=0; i<3; i++)
					// {
					// 	*deathLocationFile << cell_location[i] << "\t";
					// }
					// *deathLocationFile << node_index << "\n";
					// deathLocationFile->close();
				}
			}
		}
		

	}
}

void DensityDependantCellKillerStrip3DWithGhostNodesSimpleV2::SetNumberCellsRemovedByAnoikis(std::vector<c_vector<unsigned,2> > cellsRemoved)
{
	unsigned num_removed_by_anoikis = 0;
	
    for (unsigned i=0; i<cellsRemoved.size(); i++)
    {
    	if(cellsRemoved[i][1]==1)
    	{
    		num_removed_by_anoikis+=1;
    	}
    }

    mCellsRemovedByAnoikis += num_removed_by_anoikis;
}

unsigned DensityDependantCellKillerStrip3DWithGhostNodesSimpleV2::GetNumberCellsRemovedByAnoikis()
{
	return mCellsRemovedByAnoikis;
}

/* Data stored: time - node_index - x - y - z */
void DensityDependantCellKillerStrip3DWithGhostNodesSimpleV2::SetLocationsOfCellsRemovedByAnoikis(std::vector<c_vector<unsigned,2> > cellsRemoved)
{
	DomMeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*> (this->mpCellPopulation);
	double x_location, y_location, z_location;
	c_vector<double, 5> node_time_and_location;

	// Need to use the node indices to store the locations of where cells are removed
    for (unsigned i=0; i<cellsRemoved.size(); i++)
    {
    	if (cellsRemoved[i][1] == 1)		// This cell has been removed by anoikis
    	{
    		node_time_and_location[0] = SimulationTime::Instance()->GetTime();
    		node_time_and_location[1] = cellsRemoved[i][0];
    		
			unsigned node_index = cellsRemoved[i][0];

			CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(node_index);
			x_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[0];
			y_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[1];
			z_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[2];

			node_time_and_location[2] = x_location;
			node_time_and_location[3] = y_location;
			node_time_and_location[4] = z_location;

			mLocationsOfAnoikisCells.push_back(node_time_and_location);
    	}
    }
}

std::vector<c_vector<double,5> > DensityDependantCellKillerStrip3DWithGhostNodesSimpleV2::GetLocationsOfCellsRemovedByAnoikis()
{
	return mLocationsOfAnoikisCells;
}

void DensityDependantCellKillerStrip3DWithGhostNodesSimpleV2::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CellsRemovedByAnoikis>" << mCellsRemovedByAnoikis << "</CellsRemovedByAnoikis> \n";
	*rParamsFile << "\t\t\t<CutOffLength>" << mCutOffLength << "</CutOffLength> \n";
	*rParamsFile << "\t\t\t<DensityThreshold>" << mDensityThreshold << "</DensityThreshold> \n";
	
//    *rParamsFile << "\t\t\t<XLocationsOfAnoikisCells>" << mXLocationsOfAnoikisCells << "</XLocationsOfAnoikisCells> \n";

    // Call direct parent class
    AbstractCellKiller<3>::OutputCellKillerParameters(rParamsFile);
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(DensityDependantCellKillerStrip3DWithGhostNodesSimpleV2)
