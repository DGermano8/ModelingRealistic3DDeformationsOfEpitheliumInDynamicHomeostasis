#include "SloughingCellKiller3DWithGhostNodes.hpp"
#include "Debug.hpp"

/* Sloughing occurring for all those epithelial cells that move beyond the edges of the periodic domain (would be moved to the opposite
 * edge under complete periodic boundary conditions
 */
SloughingCellKiller3DWithGhostNodes::SloughingCellKiller3DWithGhostNodes(AbstractCellPopulation<3>* pCrypt, double cellPopulationWidth, double cellPopulationDepth)
    : AbstractCellKiller<3>(pCrypt),
      mCellPopulationWidth(cellPopulationWidth),
      mCellPopulationDepth(cellPopulationDepth),
      mCellsRemovedBySloughing(0)
{
    // Sets up output file
//	OutputFileHandler output_file_handler(mOutputDirectory + "AnoikisData/", false);
//	mAnoikisOutputFile = output_file_handler.OpenOutputFile("results.anoikis");
}

SloughingCellKiller3DWithGhostNodes::~SloughingCellKiller3DWithGhostNodes()
{
//    mAnoikisOutputFile->close();
}

double SloughingCellKiller3DWithGhostNodes::GetCellPopulationWidth() const
{
	return mCellPopulationWidth;
}

double SloughingCellKiller3DWithGhostNodes::GetCellPopulationDepth() const
{
	return mCellPopulationDepth;
}

void SloughingCellKiller3DWithGhostNodes::SetOutputDirectory(std::string outputDirectory)
{
	mOutputDirectory = outputDirectory;
}

std::string SloughingCellKiller3DWithGhostNodes::GetOutputDirectory()
{
	return mOutputDirectory;
}

/** A method to return a vector that indicates which cells should be killed by anoikis
 * and which by compression-driven apoptosis
 */
std::vector<c_vector<unsigned,2> > SloughingCellKiller3DWithGhostNodes::RemoveBySloughing()
{
	DomMeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*> (this->mpCellPopulation);

//This needs fixing
//    assert(p_tissue->GetVoronoiTessellation()!=NULL);	// This fails during archiving of a simulation as Voronoi stuff not archived yet

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

			// Examine each epithelial node to see if it should be removed by sloughing

			double x = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter)[0];
			double y = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter)[1];
			
			if ( (cell_iter->GetMutationState()->IsType<StromalCellMutationState>()==false) &&
					( (x < -0.5) || (x>mCellPopulationWidth) || (y < -0.5) || (y>mCellPopulationDepth) ) )
			{
				individual_node_information[1] = 1;
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
void SloughingCellKiller3DWithGhostNodes::CheckAndLabelCellsForApoptosisOrDeath()
{
	DomMeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*> (this->mpCellPopulation);

	//This needs fixing	
	//    assert(p_tissue->GetVoronoiTessellation()!=NULL);	// This fails during archiving of a simulation as Voronoi stuff not archived yet

    // Get the information at this timestep for each node index that says whether to remove by anoikis or sloughing
    std::vector<c_vector<unsigned,2> > cells_to_remove = this->RemoveBySloughing();

    // Keep a record of how many cells have been removed at this timestep
    this->SetNumberCellsRemovedBySloughing(cells_to_remove);
    this->SetLocationsOfCellsRemovedBySloughing(cells_to_remove);

    for (unsigned i=0; i<cells_to_remove.size(); i++)
    {
    	if (cells_to_remove[i][1] == 1)
    	{
    		// Get cell associated to this node
    		CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(cells_to_remove[i][0]);
    		p_cell->Kill();
    	}
    }
}


void SloughingCellKiller3DWithGhostNodes::SetNumberCellsRemovedBySloughing(std::vector<c_vector<unsigned,2> > cellsRemoved)
{
	unsigned num_removed_by_sloughing = 0;

    for (unsigned i=0; i<cellsRemoved.size(); i++)
    {
    	if(cellsRemoved[i][1]==1)
    	{
    		num_removed_by_sloughing+=1;
    	}
    }

    mCellsRemovedBySloughing += num_removed_by_sloughing;
}

unsigned SloughingCellKiller3DWithGhostNodes::GetNumberCellsRemovedBySloughing()
{
	return mCellsRemovedBySloughing;
}

void SloughingCellKiller3DWithGhostNodes::SetLocationsOfCellsRemovedBySloughing(std::vector<c_vector<unsigned,2> > cellsRemoved)
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

			mLocationsOfSloughedCells.push_back(node_time_and_location);
    	}
    }
}

std::vector<c_vector<double,5> > SloughingCellKiller3DWithGhostNodes::GetLocationsOfCellsRemovedBySloughing()
{
	return mLocationsOfSloughedCells;
}

void SloughingCellKiller3DWithGhostNodes::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CellsRemovedBySloughing>" << mCellsRemovedBySloughing << "</CellsRemovedBySloughing> \n";
//    *rParamsFile << "\t\t\t<LocationsOfSloughedCells>" << mLocationsOfSloughedCells << "</LocationsOfSloughedCells> \n";

    // Call direct parent class
    AbstractCellKiller<3>::OutputCellKillerParameters(rParamsFile);
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(SloughingCellKiller3DWithGhostNodes)
