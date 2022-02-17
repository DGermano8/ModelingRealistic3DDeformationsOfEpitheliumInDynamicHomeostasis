#include "UniformCellKiller3dWithGhostNodes.hpp"
#include "Debug.hpp"
#include "SimulationTime.hpp"
#include "OutputFileHandler.hpp"



/* Apoptosis for cells that are epithelial and lose contact with the basement membrane
 *
 */
UniformCellKiller3dWithGhostNodes::UniformCellKiller3dWithGhostNodes(AbstractCellPopulation<3>* pCrypt, double probabilityOfDeathInAnHour,
		double minXBoundary, double maxXBoundary, double minYBoundary, double maxYBoundary, double cellPop)
    : AbstractCellKiller<3>(pCrypt),
      mProbabilityOfDeathInAnHour(probabilityOfDeathInAnHour),
      mMinXBoundary(minXBoundary),
      mMaxXBoundary(maxXBoundary),
      mMinYBoundary(minYBoundary),
      mMaxYBoundary(maxYBoundary),
	  mMinCellPopulation(cellPop)
{
}

UniformCellKiller3dWithGhostNodes::~UniformCellKiller3dWithGhostNodes()
{
}

void UniformCellKiller3dWithGhostNodes::SetOutputDirectory(std::string outputDirectory)
{
	mOutputDirectory = outputDirectory;

	OutputFileHandler output_file_handler(mOutputDirectory+"/", false);
    out_stream deathLocationFile = output_file_handler.OpenOutputFile("uniformDeaths.dat");
    *deathLocationFile << "time \t";
    for (unsigned i=0; i<3; i++)
    {
        *deathLocationFile << "location" << i << "\t";
    }
    *deathLocationFile << "Cell ID " << "\n";
    deathLocationFile->close();
}

std::string UniformCellKiller3dWithGhostNodes::GetOutputDirectory()
{
	return mOutputDirectory;
}

/** A method to return a vector that indicates which cells should be killed by anoikis
 * and which by compression-driven apoptosis
 */
std::vector<c_vector<unsigned,2> > UniformCellKiller3dWithGhostNodes::RemoveByRandomSelection()
{
	DomMeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*> (this->mpCellPopulation);
    // assert(p_tissue->GetVoronoiTessellation()!=NULL);	// This fails during archiving of a simulation as Voronoi stuff not archived yet

    std::vector<c_vector<unsigned,2> > cells_to_remove;
    c_vector<unsigned,2> individual_node_information;	// Will store the node index and whether to remove or not (1 or 0)

	for (AbstractCellPopulation<3>::Iterator cell_iter = p_tissue->Begin();
    	 cell_iter != p_tissue->End();
    	 ++cell_iter)
	{
		unsigned node_index = p_tissue->GetNodeCorrespondingToCell(*cell_iter)->GetIndex();
		assert((!p_tissue->IsGhostNode(node_index)));

		// Initialise
		individual_node_information[0] = node_index;
		individual_node_information[1] = 0;

	    // We assume a constant time step
	    double death_prob_this_timestep = 1.0 - pow((1.0 - mProbabilityOfDeathInAnHour), SimulationTime::Instance()->GetTimeStep());
	    double x = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter)[0];
	    double y = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter)[1];

	    if ((!cell_iter->HasApoptosisBegun() && (cell_iter->GetMutationState()->IsType<StromalCellMutationState>()==false)
	    		&& (RandomNumberGenerator::Instance()->ranf() < death_prob_this_timestep) )
	    		&& ( ((x<mMinXBoundary) || (x>mMaxXBoundary)) || ((y<mMinYBoundary) || (y>mMaxYBoundary)) ) )
	    {
	    	// cell_iter->StartApoptosis();
			// PRINT_3_VARIABLES(cell_iter->GetApoptosisTime(),cell_iter->GetTimeUntilDeath(),SimulationTime::Instance()->GetTime());

	    	individual_node_information[1] = 1;
	    }

		cells_to_remove.push_back(individual_node_information);
	}

	return cells_to_remove;
}

void UniformCellKiller3dWithGhostNodes::CheckAndLabelSingleCellForApoptosis(CellPtr pCell)
{
    // We assume a constant time step
    double death_prob_this_timestep = 1.0 - pow((1.0 - mProbabilityOfDeathInAnHour), SimulationTime::Instance()->GetTimeStep());

    if (!pCell->HasApoptosisBegun() &&
        RandomNumberGenerator::Instance()->ranf() < death_prob_this_timestep)
    {
        pCell->StartApoptosis();
    }
}

/* Cell Killer that kills epithelial cells in the target zone (bounded by the specified x and y boundaries) by random selection
*/
void UniformCellKiller3dWithGhostNodes::CheckAndLabelCellsForApoptosisOrDeath()
{

	// TRACE("Entering")

	DomMeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*> (this->mpCellPopulation);
	//assert(p_tissue->GetVoronoiTessellation()!=NULL);	// This fails during archiving of a simulation as Voronoi stuff not archived yet
    double death_prob_this_timestep = 1.0 - pow((1.0 - mProbabilityOfDeathInAnHour), SimulationTime::Instance()->GetTimeStep());
	// int how_many_cells_to_kill = 0;

	// for (AbstractCellPopulation<3>::Iterator cell_iter = p_tissue->Begin();
    // 	 cell_iter != p_tissue->End();
    // 	 ++cell_iter)
	// {
	// 	unsigned node_index = p_tissue->GetNodeCorrespondingToCell(*cell_iter)->GetIndex();
	// 	assert((!p_tissue->IsGhostNode(node_index)));
	// 	CellPtr p_cell = mpCellPopulation->GetCellUsingLocationIndex(node_index);
		
	// 	// PRINT_2_VARIABLES(node_index, p_cell->GetAge());


	// 	if(p_cell->GetAge() == SimulationTime::Instance()->GetTimeStep() )
	// 	{
	// 		how_many_cells_to_kill++;
	// 	}

	// }
	// how_many_cells_to_kill = 0.5*how_many_cells_to_kill;
	// PRINT_2_VARIABLES(SimulationTime::Instance()->GetTime(),mpCellPopulation->GetNumRealCells());

	int number_of_apoptic_cells = 0;
	for (AbstractCellPopulation<3>::Iterator cell_iter = p_tissue->Begin();
    	 cell_iter != p_tissue->End();
    	 ++cell_iter)
	{
		unsigned node_index = p_tissue->GetNodeCorrespondingToCell(*cell_iter)->GetIndex();
		CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(node_index);

		if( (p_cell->GetMutationState()->IsType<WildTypeCellMutationState>()==true) && (p_cell->HasApoptosisBegun()) )
		{
			number_of_apoptic_cells++;
		}
	}

	bool can_i_kill_cells = (mpCellPopulation->GetNumRealCells() - number_of_apoptic_cells  > mMinCellPopulation);
    // Loop over this vector and kill any cells that it tells you to
    // while (how_many_cells_to_kill > 0)
	int num_cells_have_i_killed = 0;
	std::vector<unsigned> cells_which_can_die;
	
	if(can_i_kill_cells)
    {
		// need to do a random permutation here, otherwise we have some of the cells (at start of  list) having overall higher probability
		for (AbstractCellPopulation<3>::Iterator cell_iter = p_tissue->Begin();
    	 cell_iter != p_tissue->End();
    	 ++cell_iter)
		{
			unsigned node_index = p_tissue->GetNodeCorrespondingToCell(*cell_iter)->GetIndex();
			CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(node_index);
			AbstractPhaseBasedCellCycleModel* p_model = static_cast<AbstractPhaseBasedCellCycleModel*>(p_cell->GetCellCycleModel());

			if (!(p_cell->IsDead()) && !p_tissue->IsGhostNode(node_index) && 
				(p_cell->GetAge() > p_model->GetMDuration()) &&
				(p_cell->GetMutationState()->IsType<WildTypeCellMutationState>()==true) )
			{
				if(!p_cell->HasApoptosisBegun() )
				{

					double x = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter)[0];
					double y = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter)[1];

					if ( ((x<mMinXBoundary) || (x>mMaxXBoundary)) || ((y<mMinYBoundary) || (y>mMaxYBoundary)) )
					{
						
						if( !(p_cell->HasApoptosisBegun()) )
						{
							cells_which_can_die.push_back(node_index);
						}
						// cells_which_can_die.push_back(node_index);
						// if( (mpCellPopulation->GetNumRealCells() -  num_cells_have_i_killed) > mMinCellPopulation)
						// {
						// 	p_cell->StartApoptosis();
						// 	// how_many_cells_to_kill--;
						// 	num_cells_have_i_killed++;
						// 	// break;
						// }
					}

				}
			}
		}
	}

	// PRINT_3_VARIABLES(p_tissue->GetNumRealCells(),mMinCellPopulation,cells_which_can_die.size());

	std::vector<unsigned> randon_cells;
	if(can_i_kill_cells)
	{
		for(int j=0; j< (mpCellPopulation->GetNumRealCells() - number_of_apoptic_cells - mMinCellPopulation); j++)
		{
			int cell_to_kill = int ((cells_which_can_die.size() )*(RandomNumberGenerator::Instance()->ranf()));

			// PRINT_3_VARIABLES(cell_to_kill,mpCellPopulation->GetNumRealCells(),cells_which_can_die.size());

			CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(cells_which_can_die[cell_to_kill]);

			if( !(p_cell->HasApoptosisBegun()) )
			{
				p_cell->StartApoptosis();
			}
			// TRACE("Killed");
		}
	}

	// The actual killer
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

					SimulationTime* p_time = SimulationTime::Instance();
            		c_vector<double, 3> cell_location = p_tissue->GetNode(node_index)->rGetLocation();


					OutputFileHandler output_file_handler(mOutputDirectory+"/", false);
					out_stream deathLocationFile = output_file_handler.OpenOutputFile("uniformDeaths.dat", std::ios::app);

					*deathLocationFile << p_time->GetTime() << "\t";
					for (unsigned i=0; i<3; i++)
					{
						*deathLocationFile << cell_location[i] << "\t";
					}
					*deathLocationFile << node_index << "\n";
					deathLocationFile->close();
				}
			}
		}
		
	}
	// TRACE("Leaving")

}

void UniformCellKiller3dWithGhostNodes::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ProbabilityOfDeathInAnHour>" << mProbabilityOfDeathInAnHour << "</ProbabilityOfDeathInAnHour> \n";
    *rParamsFile << "\t\t\t<MinXBoundary>" << mMinXBoundary << "</MinXBoundary> \n";
    *rParamsFile << "\t\t\t<MaxXBoundary>" << mMaxXBoundary << "</MaxXBoundary> \n";
    *rParamsFile << "\t\t\t<MinYBoundary>" << mMinYBoundary << "</MinYBoundary> \n";
    *rParamsFile << "\t\t\t<MaxYBoundary>" << mMaxYBoundary << "</MaxYBoundary> \n";
    *rParamsFile << "\t\t\t<CellPopulation>" << mMinCellPopulation << "</CellPopulation> \n";

	

    // Call direct parent class
    AbstractCellKiller<3>::OutputCellKillerParameters(rParamsFile);
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(UniformCellKiller3dWithGhostNodes)
