#include "PeriodicBoxBoundaryCondition3dGhosts.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "StromalCellMutationState.hpp"
#include "Debug.hpp"

PeriodicBoxBoundaryCondition3dGhosts::PeriodicBoxBoundaryCondition3dGhosts(AbstractCellPopulation<3>* pCellPopulation)
    : AbstractCellPopulationBoundaryCondition<3>(pCellPopulation),
    mCellPopulationWidth(DOUBLE_UNSET),
    mCellPopulationDepth(DOUBLE_UNSET),
    mMaxHeightForPinnedCells(DOUBLE_UNSET)
{
}

PeriodicBoxBoundaryCondition3dGhosts::~PeriodicBoxBoundaryCondition3dGhosts()
{
} 


//void PeriodicBoxBoundaryCondition3d::ImposeBoundaryCondition(const std::vector< c_vector<double, 3> >& rOldLocations)
void PeriodicBoxBoundaryCondition3dGhosts::ImposeBoundaryCondition(const std::map<Node<3>*, c_vector<double, 3> >& rOldLocations)
{	
	// TRACE("Imposing BCs");
	// for (AbstractCellPopulation<3>::Iterator cell_iter = this->mpCellPopulation->Begin();
	// 	 cell_iter != this->mpCellPopulation->End();
	// 	 ++cell_iter)
	for(unsigned node_index=0; node_index<this->mpCellPopulation->GetNumNodes(); node_index++)
	{
		
		// assert(this->mpCellPopulation->IsCellAssociatedWithADeletedLocation(*cell_iter) == false);

		// Get index of node associated with cell
        // unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);

        // Get pointer to this node
        Node<3>* p_node = this->mpCellPopulation->GetNode(node_index);

		// Pin the bottom cells (only stromal)
		// if ( (this->mpCellPopulation->GetCellUsingLocationIndex(node_index)->GetMutationState()->IsType<StromalCellMutationState>())
		// 		&& (p_node->rGetLocation()[2] < mMaxHeightForPinnedCells) )
		// {
		// 	// Get old node location
		// 	typename std::map<Node<3>*, c_vector<double, 3> >::const_iterator it = rOldLocations.find(p_node);
    	// 	c_vector<double, 3> old_node_location = it->second;
		// 	// c_vector<double, 3> old_node_location = rOldLocations[node_index];

		// 	// Return node to old location
		// 	p_node->rGetModifiableLocation()[0] = old_node_location[0];
		// 	p_node->rGetModifiableLocation()[1] = old_node_location[1];
		// 	// p_node->rGetModifiableLocation()[2] = 0.0;
		// 	p_node->rGetModifiableLocation()[2] = mMaxHeightForPinnedCells;
		// }

		// Want to enforce a periodic boundary box, such that cells which move outside the boundaries of the box get
		// moved to the opposite side (as in the cylindrical mesh)
		if (p_node->rGetLocation()[0] < 0.0)
		{
			/* Move to the opposite edge */
			 p_node->rGetModifiableLocation()[0] = p_node->rGetLocation()[0] + mCellPopulationWidth;
		}
		else if (p_node->rGetLocation()[0] > mCellPopulationWidth)
		{
			/* Move to the opposite edge */
			p_node->rGetModifiableLocation()[0] = p_node->rGetLocation()[0] - mCellPopulationWidth;
		}
		if (p_node->rGetLocation()[1] < 0.0)
		{
			/* Move to the opposite edge */
			 p_node->rGetModifiableLocation()[1] = p_node->rGetLocation()[1] + mCellPopulationDepth;
		}
		else if (p_node->rGetLocation()[1] > mCellPopulationDepth)
		{
			/* Move to the opposite edge */
			p_node->rGetModifiableLocation()[1] = p_node->rGetLocation()[1] - mCellPopulationDepth;
		}
	}
}

bool PeriodicBoxBoundaryCondition3dGhosts::VerifyBoundaryCondition()
{
	bool boundary_condition_satisfied = true;

	/*
	 * Here we verify that the boundary condition is still satisfied by simply
	 * checking that no cells lies below the Z=0 boundary.
	 */
    // for (AbstractCellPopulation<3>::Iterator cell_iter = this->mpCellPopulation->Begin();
    //      cell_iter != this->mpCellPopulation->End();
    //      ++cell_iter)
	for(unsigned node_index=0; node_index<this->mpCellPopulation->GetNumNodes(); node_index++)
    {
        // Get index of node associated with cell
        // unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);

        // Get pointer to this node
        Node<3>* p_node = this->mpCellPopulation->GetNode(node_index);

        // If this node lies below the z=0 boundary (NO! - Dom), or outside of the periodic boundaries, break and return false
        // If this node lies below the z=0 boundary, or outside of the periodic boundaries, break and return false
        // if ( (p_node->rGetLocation()[2] < 0.0) || (p_node->rGetLocation()[0] < 0.0) || (p_node->rGetLocation()[0] > mCellPopulationWidth)
        // 		|| (p_node->rGetLocation()[1] < 0.0) || (p_node->rGetLocation()[1] > mCellPopulationDepth) )
        if ( (p_node->rGetLocation()[0] < 0.0) || (p_node->rGetLocation()[0] > mCellPopulationWidth)
          || (p_node->rGetLocation()[1] < 0.0) || (p_node->rGetLocation()[1] > mCellPopulationDepth) )
        {
        	PRINT_2_VARIABLES(node_index, SimulationTime::Instance()->GetTime());
        	PRINT_VECTOR(p_node->rGetLocation());
        	boundary_condition_satisfied = false;
        	break;
        }
    }

    return boundary_condition_satisfied;
}

void PeriodicBoxBoundaryCondition3dGhosts::SetMaxHeightForPinnedCells(double maxHeightForPinnedCells)
{
	mMaxHeightForPinnedCells = maxHeightForPinnedCells;
}

double PeriodicBoxBoundaryCondition3dGhosts::GetMaxHeightForPinnedCells()
{
	return mMaxHeightForPinnedCells;
}

void PeriodicBoxBoundaryCondition3dGhosts::SetCellPopulationWidth(double cellPopulationWidth)
{
	mCellPopulationWidth = cellPopulationWidth;
}

double PeriodicBoxBoundaryCondition3dGhosts::GetCellPopulationWidth()
{
	return mCellPopulationWidth;
}

void PeriodicBoxBoundaryCondition3dGhosts::SetCellPopulationDepth(double cellPopulationDepth)
{
	mCellPopulationDepth = cellPopulationDepth;
}

double PeriodicBoxBoundaryCondition3dGhosts::GetCellPopulationDepth()
{
	return mCellPopulationDepth;
}


void PeriodicBoxBoundaryCondition3dGhosts::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<CellPopulationDepth>" << mCellPopulationDepth << "</CellPopulationDepth>\n";
    *rParamsFile << "\t\t<CellPopulationWidth>" << mCellPopulationWidth << "</CellPopulationWidth>\n";
    *rParamsFile << "\t\t<MaxHeightForPinnedCells>" << mMaxHeightForPinnedCells << "</MaxHeightForPinnedCells>\n";
    
    // Call method on direct parent class
	AbstractCellPopulationBoundaryCondition<3>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PeriodicBoxBoundaryCondition3dGhosts)
