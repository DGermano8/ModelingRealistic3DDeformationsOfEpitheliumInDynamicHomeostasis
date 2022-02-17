#include "PeriodicStromalBoxBoundaryCondition3d.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "StromalCellMutationState.hpp"

PeriodicStromalBoxBoundaryCondition3d::PeriodicStromalBoxBoundaryCondition3d(AbstractCellPopulation<3>* pCellPopulation)
    : AbstractCellPopulationBoundaryCondition<3>(pCellPopulation),
    mCellPopulationWidth(DOUBLE_UNSET),
    mCellPopulationDepth(DOUBLE_UNSET),
    mMaxHeightForPinnedCells(DOUBLE_UNSET)
{
}

PeriodicStromalBoxBoundaryCondition3d::~PeriodicStromalBoxBoundaryCondition3d()
{
}
 //void PeriodicStromalBoxBoundaryCondition3d::ImposeBoundaryCondition(const std::vector< c_vector<double, 3> >& rOldLocations)
void PeriodicStromalBoxBoundaryCondition3d::ImposeBoundaryCondition(const std::map<Node<3>*, c_vector<double, 3> >& rOldLocations)
{	
	for (AbstractCellPopulation<3>::Iterator cell_iter = this->mpCellPopulation->Begin();
		 cell_iter != this->mpCellPopulation->End();
		 ++cell_iter)
	{
		assert(this->mpCellPopulation->IsCellAssociatedWithADeletedLocation(*cell_iter) == false);

		// Get index of node associated with cell
        unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);

        // Get pointer to this node
        Node<3>* p_node = this->mpCellPopulation->GetNode(node_index);
		
		if (this->mpCellPopulation->GetCellUsingLocationIndex(node_index)->GetMutationState()->IsType<StromalCellMutationState>())
		{
			// Pin the bottom cells 
			if (p_node->rGetLocation()[2] < mMaxHeightForPinnedCells) 
			{
				// Get old node location
				//c_vector<double, 3> old_node_location = rOldLocations[node_index];
				typename std::map<Node<3>*, c_vector<double, 3> >::const_iterator it = rOldLocations.find(p_node);
    			c_vector<double, 3> old_node_location = it->second;
	
				// Return node to old location
				p_node->rGetModifiableLocation()[0] = old_node_location[0];
				p_node->rGetModifiableLocation()[1] = old_node_location[1];
				// p_node->rGetModifiableLocation()[2] = 0.0;
				p_node->rGetModifiableLocation()[2] = mMaxHeightForPinnedCells;
			}
	
			// Want to enforce a periodic boundary box on only the stromal cells, such that cells which move outside the 
			// boundaries of the box get moved to the opposite side (as in the cylindrical mesh)
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
}

bool PeriodicStromalBoxBoundaryCondition3d::VerifyBoundaryCondition()
{
	bool boundary_condition_satisfied = true;

	/*
	 * Here we verify that the boundary condition is still satisfied by simply
	 * checking that no cells lies below the Z=0 boundary.
	 */
    for (AbstractCellPopulation<3>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        // Get index of node associated with cell
        unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);

        // Get pointer to this node
        Node<3>* p_node = this->mpCellPopulation->GetNode(node_index);

        // If this node lies below the z=0 boundary(NO! - Dom), or outside of the periodic boundaries, break and return false
		// if ( (cell_iter->GetMutationState()->IsType<StromalCellMutationState>()) && 
        // 		((p_node->rGetLocation()[2] < 0.0) || (p_node->rGetLocation()[0] < 0.0) || (p_node->rGetLocation()[0] > mCellPopulationWidth)
        // 		|| (p_node->rGetLocation()[1] < 0.0) || (p_node->rGetLocation()[1] > mCellPopulationDepth)) )

        // If this node is outside of the periodic boundaries, break and return false
        if ( (cell_iter->GetMutationState()->IsType<StromalCellMutationState>()) && 
        		  ((p_node->rGetLocation()[0] < 0.0) || (p_node->rGetLocation()[0] > mCellPopulationWidth)
        		|| (p_node->rGetLocation()[1] < 0.0) || (p_node->rGetLocation()[1] > mCellPopulationDepth)) )
        {
        	boundary_condition_satisfied = false;
        	break;
        }
    }

    return boundary_condition_satisfied;
}

void PeriodicStromalBoxBoundaryCondition3d::SetMaxHeightForPinnedCells(double maxHeightForPinnedCells)
{
	mMaxHeightForPinnedCells = maxHeightForPinnedCells;
}

double PeriodicStromalBoxBoundaryCondition3d::GetMaxHeightForPinnedCells()
{
	return mMaxHeightForPinnedCells;
}

void PeriodicStromalBoxBoundaryCondition3d::SetCellPopulationWidth(double cellPopulationWidth)
{
	mCellPopulationWidth = cellPopulationWidth;
}

double PeriodicStromalBoxBoundaryCondition3d::GetCellPopulationWidth()
{
	return mCellPopulationWidth;
}

void PeriodicStromalBoxBoundaryCondition3d::SetCellPopulationDepth(double cellPopulationDepth)
{
	mCellPopulationDepth = cellPopulationDepth;
}

double PeriodicStromalBoxBoundaryCondition3d::GetCellPopulationDepth()
{
	return mCellPopulationDepth;
}


void PeriodicStromalBoxBoundaryCondition3d::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<CellPopulationDepth>" << mCellPopulationDepth << "</CellPopulationDepth>\n";
    *rParamsFile << "\t\t<CellPopulationWidth>" << mCellPopulationWidth << "</CellPopulationWidth>\n";
    *rParamsFile << "\t\t<MaxHeightForPinnedCells>" << mMaxHeightForPinnedCells << "</MaxHeightForPinnedCells>\n";
    
    // Call method on direct parent class
	AbstractCellPopulationBoundaryCondition<3>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PeriodicStromalBoxBoundaryCondition3d)
