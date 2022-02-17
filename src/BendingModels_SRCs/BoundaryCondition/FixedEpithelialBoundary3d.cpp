#include "FixedEpithelialBoundary3d.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "WildTypeCellMutationState.hpp"

FixedEpithelialBoundary3d::FixedEpithelialBoundary3d(AbstractCellPopulation<3>* pCellPopulation)
    : AbstractCellPopulationBoundaryCondition<3>(pCellPopulation),
    mCellPopulationWidth(DOUBLE_UNSET),
    mCellPopulationDepth(DOUBLE_UNSET),
    mHeightForPinnedCells(DOUBLE_UNSET)
{
}

FixedEpithelialBoundary3d::~FixedEpithelialBoundary3d()
{
}
 //void FixedEpithelialBoundary3d::ImposeBoundaryCondition(const std::vector< c_vector<double, 3> >& rOldLocations)
void FixedEpithelialBoundary3d::ImposeBoundaryCondition(const std::map<Node<3>*, c_vector<double, 3> >& rOldLocations)
{	
	double eps = 0.0000001;

	for (AbstractCellPopulation<3>::Iterator cell_iter = this->mpCellPopulation->Begin();
		 cell_iter != this->mpCellPopulation->End();
		 ++cell_iter)
	{
		assert(this->mpCellPopulation->IsCellAssociatedWithADeletedLocation(*cell_iter) == false);

		// Get index of node associated with cell
        unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);

        // Get pointer to this node
        Node<3>* p_node = this->mpCellPopulation->GetNode(node_index);
		
		if (this->mpCellPopulation->GetCellUsingLocationIndex(node_index)->GetMutationState()->IsType<WildTypeCellMutationState>())
		{
			// Pin the epithelial cells 
			if ( (p_node->rGetLocation()[2] < mHeightForPinnedCells - eps) || (p_node->rGetLocation()[2] > mHeightForPinnedCells + eps) )
			{
				// Get old node location
				//c_vector<double, 3> old_node_location = rOldLocations[node_index];
				// typename std::map<Node<3>*, c_vector<double, 3> >::const_iterator it = rOldLocations.find(p_node);
    			// c_vector<double, 3> old_node_location = it->second;
	
				// Return node to old location
				p_node->rGetModifiableLocation()[0] = p_node->rGetLocation()[0];
				p_node->rGetModifiableLocation()[1] = p_node->rGetLocation()[1];
				// p_node->rGetModifiableLocation()[2] = 0.0;
				p_node->rGetModifiableLocation()[2] = mHeightForPinnedCells;
			}
		}
	}
}

bool FixedEpithelialBoundary3d::VerifyBoundaryCondition()
{
	bool boundary_condition_satisfied = true;
	double eps = 0.000001;

	/*
	 * Here we verify that the boundary condition is still satisfied by simply
	 * checking that no cells lies below the Z=0 boundary.
	 */
    for (AbstractCellPopulation<3>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
		assert(this->mpCellPopulation->IsCellAssociatedWithADeletedLocation(*cell_iter) == false);
		
        // Get index of node associated with cell
        unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);

        // Get pointer to this node
        Node<3>* p_node = this->mpCellPopulation->GetNode(node_index);

        // If this node lies below the z=0 boundary(NO! - Dom), or outside of the periodic boundaries, break and return false
		// if ( (cell_iter->GetMutationState()->IsType<StromalCellMutationState>()) && 
        // 		((p_node->rGetLocation()[2] < 0.0) || (p_node->rGetLocation()[0] < 0.0) || (p_node->rGetLocation()[0] > mCellPopulationWidth)
        // 		|| (p_node->rGetLocation()[1] < 0.0) || (p_node->rGetLocation()[1] > mCellPopulationDepth)) )

        // If this node is outside of the periodic boundaries, break and return false
        if ( (cell_iter->GetMutationState()->IsType<WildTypeCellMutationState>()) &&
			(p_node->rGetLocation()[2] < mHeightForPinnedCells - eps) || (p_node->rGetLocation()[2] > mHeightForPinnedCells + eps) )
        {
        	boundary_condition_satisfied = false;
        	break;
        }
    }

    return boundary_condition_satisfied;
}

void FixedEpithelialBoundary3d::SetHeightForPinnedCells(double heightForPinnedCells)
{
	mHeightForPinnedCells = heightForPinnedCells;
}

double FixedEpithelialBoundary3d::GetHeightForPinnedCells()
{
	return mHeightForPinnedCells;
}

void FixedEpithelialBoundary3d::SetCellPopulationWidth(double cellPopulationWidth)
{
	mCellPopulationWidth = cellPopulationWidth;
}

double FixedEpithelialBoundary3d::GetCellPopulationWidth()
{
	return mCellPopulationWidth;
}

void FixedEpithelialBoundary3d::SetCellPopulationDepth(double cellPopulationDepth)
{
	mCellPopulationDepth = cellPopulationDepth;
}

double FixedEpithelialBoundary3d::GetCellPopulationDepth()
{
	return mCellPopulationDepth;
}


void FixedEpithelialBoundary3d::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<CellPopulationDepth>" << mCellPopulationDepth << "</CellPopulationDepth>\n";
    *rParamsFile << "\t\t<CellPopulationWidth>" << mCellPopulationWidth << "</CellPopulationWidth>\n";
    *rParamsFile << "\t\t<MaxHeightForPinnedCells>" << mHeightForPinnedCells << "</MaxHeightForPinnedCells>\n";
    
    // Call method on direct parent class
	AbstractCellPopulationBoundaryCondition<3>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(FixedEpithelialBoundary3d)
