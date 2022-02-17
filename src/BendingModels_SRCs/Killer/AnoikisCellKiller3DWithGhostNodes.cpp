#include "AnoikisCellKiller3DWithGhostNodes.hpp"
#include "OutputFileHandler.hpp"
#include "Debug.hpp"

/* Apoptosis for cells that are epithelial and lose contact with the basement membrane
 *
 */
AnoikisCellKiller3DWithGhostNodes::AnoikisCellKiller3DWithGhostNodes(AbstractCellPopulation<3>* pCrypt,
																	double cut_off,
																	double cellPopulationWidth,
																	double cellPopulationDepth)
    : AbstractCellKiller<3>(pCrypt),
    mCellsRemovedByAnoikis(0),
	mCutOffLength(cut_off),
	mCellPopulationWidth(cellPopulationWidth),
    mCellPopulationDepth(cellPopulationDepth)
{
    // Sets up output file
//	OutputFileHandler output_file_handler(mOutputDirectory + "AnoikisData/", false);
//	mAnoikisOutputFile = output_file_handler.OpenOutputFile("results.anoikis");
	
	// TRACE("Start of anoikus");
	// TRACE(GetOutputDirectory());
	// OutputFileHandler output_file_handler(GetOutputDirectory()+"/", false);
	// TRACE("A");
    // out_stream deathLocationFile = output_file_handler.OpenOutputFile("anoikisDeaths.dat");
	// TRACE("B");
    // *deathLocationFile << "time \t";
	// TRACE("C");
    // for (unsigned i=0; i<3; i++)
    // {
    //     *deathLocationFile << "location" << i << "\t";
    // }
	// TRACE("D");
    // *deathLocationFile << "Cell ID " << "\n";
	// TRACE("E");
    // deathLocationFile->close();
	// TRACE("closed file");
}

AnoikisCellKiller3DWithGhostNodes::~AnoikisCellKiller3DWithGhostNodes()
{
//    mAnoikisOutputFile->close();
}

void AnoikisCellKiller3DWithGhostNodes::SetOutputDirectory(std::string outputDirectory)
{
	mOutputDirectory = outputDirectory;
	
	OutputFileHandler output_file_handler(mOutputDirectory+"/", false);
    out_stream deathLocationFile = output_file_handler.OpenOutputFile("anoikisDeaths.dat");
    *deathLocationFile << "time \t";
    for (unsigned i=0; i<3; i++)
    {
        *deathLocationFile << "location" << i << "\t";
    }
    *deathLocationFile << "Cell ID " << "\n";
    deathLocationFile->close();
}

std::string AnoikisCellKiller3DWithGhostNodes::GetOutputDirectory()
{
	return mOutputDirectory;
}

/*
 * Method to get the neighbouring nodes (excluding ghost nodes) of a particular node
 * Can then be used to identify the type of cells that surround a particular cell.
 */
std::set<unsigned> AnoikisCellKiller3DWithGhostNodes::GetNeighbouringNodeIndices(unsigned nodeIndex)
{
	// Create a set of neighbouring node indices
	std::set<unsigned> neighbouring_node_indices;

	// Need access to the mesh but can't get to it because the cell killer only owns a
	// pointer to an AbstractCellPopulation
    DomMeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*> (this->mpCellPopulation);

	// Find the indices of the elements owned by this node
    std::set<unsigned> containing_elem_indices = p_tissue->rGetMesh().GetNode(nodeIndex)->rGetContainingElementIndices();

	CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(nodeIndex);
	double x_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[0];
	double y_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[1];
	double z_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[2];

	//double cut_off = 1.5;

    // Iterate over these elements
    for (std::set<unsigned>::iterator elem_iter=containing_elem_indices.begin();
         elem_iter != containing_elem_indices.end();
	     ++elem_iter)
    {
	    // Get all the nodes contained in this element
	    unsigned neighbour_global_index;

	    for (unsigned local_index=0; local_index<4; local_index++)
	    {
	    	neighbour_global_index = p_tissue->rGetMesh().GetElement(*elem_iter)->GetNodeGlobalIndex(local_index);
	    	// Don't want to include the original node or ghost nodes
			if(!p_tissue->IsGhostNode(neighbour_global_index))
			{
				CellPtr p_cell_n = p_tissue->GetCellUsingLocationIndex(neighbour_global_index);
				double x_location_n = this->mpCellPopulation->GetLocationOfCellCentre(p_cell_n)[0];
				double y_location_n = this->mpCellPopulation->GetLocationOfCellCentre(p_cell_n)[1];
				double z_location_n = this->mpCellPopulation->GetLocationOfCellCentre(p_cell_n)[2];

				bool less_than_cutt_off = pow((x_location-x_location_n),2)+pow((y_location-y_location_n),2)+pow((z_location-z_location_n),2) < pow(mCutOffLength,2);

				if( (neighbour_global_index != nodeIndex) && less_than_cutt_off)
				{
					neighbouring_node_indices.insert(neighbour_global_index);
				}
			}
			else if(p_tissue->IsGhostNode(neighbour_global_index))
			{
				//const c_vector<double, 3>& node_location = this->GetNode(neighbour_global_index)->rGetLocation();
				//bool less_than_cutt_off = sqrt(pow((x_location-node_location[0]),2)+pow((y_location-node_location[1]),2)+pow((z_location-node_location[2]),2)) < mCutOffLength;

				if( (neighbour_global_index != nodeIndex))// && less_than_cutt_off)
				{
					neighbouring_node_indices.insert(neighbour_global_index);
				}
			}
	    }
    }
    return neighbouring_node_indices;		// This will contain repeats, but it doesn't matter
}

/** Method to determine if an epithelial cell has lost all contacts with the stromal cells below
 * TRUE if cell has popped up
 * FALSE if cell remains in the monolayer
 */
bool AnoikisCellKiller3DWithGhostNodes::HasCellPoppedUp(unsigned nodeIndex)
{
	DomMeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*> (this->mpCellPopulation);

	bool has_cell_popped_up = false;	// Initialising

   	std::set<unsigned> neighbours = GetNeighbouringNodeIndices(nodeIndex);

   	unsigned num_stromal_neighbours = 0;
	unsigned num_ghost_neighbours = 0;
	unsigned num_epithelial_neighbours = 0;

   	// Iterate over the neighbouring cells to check the number of stromal cell neighbours
   	for(std::set<unsigned>::iterator neighbour_iter=neighbours.begin();
   							neighbour_iter != neighbours.end();
   							++neighbour_iter)
	{
		// Check how many stromal cells its conected to
		if ( (!p_tissue->IsGhostNode(*neighbour_iter))
				&& (p_tissue->GetCellUsingLocationIndex(*neighbour_iter)->GetMutationState()->IsType<StromalCellMutationState>()==true) )
   		{
			num_stromal_neighbours += 1;
		}
		else if ( (!p_tissue->IsGhostNode(*neighbour_iter))
				&& (p_tissue->GetCellUsingLocationIndex(*neighbour_iter)->GetMutationState()->IsType<WildTypeCellMutationState>()==true) )
   		{
			num_epithelial_neighbours += 1;
		}
		// Check how many ghost nodes its conected to
		else if (p_tissue->IsGhostNode(*neighbour_iter))
		{
			num_ghost_neighbours += 1;
		}
   	}
	// PRINT_2_VARIABLES(num_stromal_neighbours,num_ghost_neighbours);
	// Ensure that the cell is connected to atleast 1 stromal cell and atleast 1 ghost node (not needed, only remove those not connected to stroma)
   	if(num_stromal_neighbours < 1 || num_epithelial_neighbours < 1 ) //|| num_ghost_neighbours < 1)
   	{
		// PRINT_2_VARIABLES(num_stromal_neighbours,num_ghost_neighbours);
		// PRINT_4_VARIABLES(SimulationTime::Instance()->GetTime(),nodeIndex, neighbours.size(), num_stromal_neighbours);
   		has_cell_popped_up = true;
		// PRINT_3_VARIABLES(SimulationTime::Instance()->GetTime(), num_stromal_neighbours,num_ghost_neighbours);
   	}
	//PRINT_VARIABLE(num_stromal_neighbours);
   
	return has_cell_popped_up;
}

/** A method to return a vector that indicates which cells should be killed by anoikis
 */
std::vector<c_vector<unsigned,2> > AnoikisCellKiller3DWithGhostNodes::RemoveByAnoikis()
{
	DomMeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*> (this->mpCellPopulation);
	double domain_tollerance = 0.1;
	
	//This needs fixing
	//   assert(p_tissue->GetVoronoiTessellation()!=NULL);	// This fails during archiving of a simulation as Voronoi stuff not archived yet

    std::vector<c_vector<unsigned,2> > cells_to_remove;
    c_vector<unsigned,2> individual_node_information;	// Will store the node index and whether to remove or not (1 or 0)

	for (AbstractCellPopulation<3>::Iterator cell_iter = p_tissue->Begin();
    	 cell_iter != p_tissue->End();
    	 ++cell_iter)
	{
		unsigned node_index = p_tissue->GetNodeCorrespondingToCell(*cell_iter)->GetIndex();
		//assert((!p_tissue->IsGhostNode(node_index)));
		if(!p_tissue->IsGhostNode(node_index))
		{
			// Initialise
			individual_node_information[0] = node_index;
			individual_node_information[1] = 0;

			// Examine each epithelial node to see if it should be removed by anoikis
			if (cell_iter->GetMutationState()->IsType<StromalCellMutationState>()==false)
			{
				assert(cell_iter->GetMutationState()->IsType<WildTypeCellMutationState>()==true);
				
				// Determining whether to remove this cell by anoikis
				if((!cell_iter->IsDead()) && (this->HasCellPoppedUp(node_index)==true) )
				{

					CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(node_index);

					// c_vector<double, 3> node_A_location = p_tissue->GetNode(node_index)->rGetLocation();

					double x_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[0];
					double y_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[1];
					double z_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[2];
					
					// Make sure cell is well inside tissue and not just an edge
					if (x_location >= domain_tollerance && x_location <= mCellPopulationWidth - domain_tollerance)
					{
						if (y_location >= domain_tollerance && y_location <= mCellPopulationDepth - domain_tollerance)
						{
							individual_node_information[1] = 1;
						}
					}
					
				}
			}

			cells_to_remove.push_back(individual_node_information);
		}
	}

	return cells_to_remove;
}

/* Cell Killer that kills transit cells that move beyond the walls of the box
 * and also any transit cells that pop upwards and become detached from the stromal
 * cells
*/
void AnoikisCellKiller3DWithGhostNodes::CheckAndLabelCellsForApoptosisOrDeath()
{
	// TRACE("Anoikis In")
	DomMeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*> (this->mpCellPopulation);

	//This needs fixing
	//    assert(p_tissue->GetVoronoiTessellation()!=NULL);	// This fails during archiving of a simulation as Voronoi stuff not archived yet

    // Get the information at this timestep for each node index that says whether to remove by anoikis 
    std::vector<c_vector<unsigned,2> > cells_to_remove = this->RemoveByAnoikis();

    // Keep a record of how many cells have been removed at this timestep
    this->SetNumberCellsRemovedByAnoikis(cells_to_remove);
    this->SetLocationsOfCellsRemovedByAnoikis(cells_to_remove);

    // Need to avoid trying to kill any cells twice (i.e. both by anoikis or random apoptosis)
    // Loop over these vectors individually and kill any cells that they tell you to

    for (unsigned i=0; i<cells_to_remove.size(); i++)
    {
    	if (cells_to_remove[i][1] == 1)
    	{
			// PRINT_VARIABLE("cell removed");
    		// Get cell associated to this node
			
    		CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(cells_to_remove[i][0]);
			p_cell->Kill();

			SimulationTime* p_time = SimulationTime::Instance();
            c_vector<double, 3> cell_location = p_tissue->GetNode(cells_to_remove[i][0])->rGetLocation();


			OutputFileHandler output_file_handler(mOutputDirectory+"/", false);
            out_stream deathLocationFile = output_file_handler.OpenOutputFile("anoikisDeaths.dat", std::ios::app);

            *deathLocationFile << p_time->GetTime() << "\t";
            for (unsigned i=0; i<3; i++)
            {
                *deathLocationFile << cell_location[i] << "\t";
            }
            *deathLocationFile << cells_to_remove[i][0] << "\n";
            deathLocationFile->close();

			// TRACE("Cell removed by Anoikis");
			
			// c_vector<double, 3> node_A_location = p_tissue->GetNode(cells_to_remove[i][0])->rGetLocation();

			// double x_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[0];
			// double y_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[1];
			// double z_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[2];

			// PRINT_3_VARIABLES(x_location,y_location,z_location);
    		
			// // Make sure cell is well inside tissue and not just an edge
			// if (x_location >= 1.0 && x_location <= mCellPopulationWidth - 1.0)
			// {
			// 	if (y_location >= 1.0 && y_location <= mCellPopulationDepth - 1.0)
			// 	{
			// 		TRACE("Yup, inside");
			// 	}
			// }
			
			
    	}
    }
	// TRACE("Anoikis Out")
}


void AnoikisCellKiller3DWithGhostNodes::SetNumberCellsRemovedByAnoikis(std::vector<c_vector<unsigned,2> > cellsRemoved)
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

unsigned AnoikisCellKiller3DWithGhostNodes::GetNumberCellsRemovedByAnoikis()
{
	return mCellsRemovedByAnoikis;
}

/* Data stored: time - node_index - x - y - z */
void AnoikisCellKiller3DWithGhostNodes::SetLocationsOfCellsRemovedByAnoikis(std::vector<c_vector<unsigned,2> > cellsRemoved)
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

std::vector<c_vector<double,5> > AnoikisCellKiller3DWithGhostNodes::GetLocationsOfCellsRemovedByAnoikis()
{
	return mLocationsOfAnoikisCells;
}

void AnoikisCellKiller3DWithGhostNodes::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CellsRemovedByAnoikis>" << mCellsRemovedByAnoikis << "</CellsRemovedByAnoikis> \n";
	*rParamsFile << "\t\t\t<CutOffLength>" << mCutOffLength << "</CutOffLength> \n";
	
//    *rParamsFile << "\t\t\t<XLocationsOfAnoikisCells>" << mXLocationsOfAnoikisCells << "</XLocationsOfAnoikisCells> \n";

    // Call direct parent class
    AbstractCellKiller<3>::OutputCellKillerParameters(rParamsFile);
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(AnoikisCellKiller3DWithGhostNodes)
