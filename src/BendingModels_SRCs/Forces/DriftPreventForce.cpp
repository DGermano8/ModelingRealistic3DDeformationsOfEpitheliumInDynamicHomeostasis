#include "DriftPreventForce.hpp"
// #include "Debug.hpp"

DriftPreventForce::DriftPreventForce()
    : AbstractForce<3>(),
	  mMovementParameter(0.001),
      mTissueMiddle(0.0)
{
}

DriftPreventForce::~DriftPreventForce()
{
}

void DriftPreventForce::SetMovementParameter(double movementParameter)
{
    assert(movementParameter > 0.0);
    mMovementParameter = movementParameter;
}

double DriftPreventForce::GetMovementParameter()
{
    return mMovementParameter;
}

void DriftPreventForce::SetTissueMiddle(double tissue_middle)
{
	mTissueMiddle = tissue_middle;
}

double DriftPreventForce::GetTissueMiddle()
{
	return mTissueMiddle;
}

void DriftPreventForce::AddForceContribution(AbstractCellPopulation<3>& rCellPopulation)
{
    // TRACE("yo");
    double dt = mMovementParameter;

   	double tissue_max_height = -10.0;
	double tissue_min_height = 10.0;
	double average_height = 0.0;

    double tissue_middle = GetTissueMiddle();
    // switch (DIM)
    // {
    //     case 1:
    //     {
    //         break;
    //     }
    //     case 2:
    //     {
    //         break;
    //     }
    //     case 3:
    //     {
            for (typename AbstractCellPopulation<3>::Iterator cell_iter = rCellPopulation.Begin();
                cell_iter != rCellPopulation.End();
                ++cell_iter)
            {
                // First, create and store a copy of this real node and cell
                unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                c_vector<double, 3> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

                if (real_node_location[2] > tissue_max_height)
                {
                    tissue_max_height = real_node_location[2];
                }
                if (real_node_location[2] < tissue_min_height)
                {
                    tissue_min_height = real_node_location[2];
                }
                average_height = average_height + real_node_location[2];

            }
            average_height = average_height/(rCellPopulation.GetNumRealCells());

            // TRACE("v1");

            // Iterate over the nodes
        
            for (typename AbstractCellPopulation<3>::Iterator cell_iter = rCellPopulation.Begin();
            cell_iter != rCellPopulation.End();
            ++cell_iter)
            {
                unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

                c_vector<double, 3> force_due_to_drift;
                force_due_to_drift[0] = 0.0; force_due_to_drift[1] = 0.0;

                force_due_to_drift[2] = -1.0*(average_height - tissue_middle)/(dt);

                rCellPopulation.GetNode(real_node_index)->AddAppliedForceContribution(force_due_to_drift);
                // node_iter->AddAppliedForceContribution(force_due_to_drift);
                // PRINT_VECTOR(force_due_to_drift);
            }
            // TRACE("v2");
    //     }
        
    // }
}

void DriftPreventForce::OutputForceParameters(out_stream& rParamsFile)
{
    // *rParamsFile << "\t\t\t<MovementParameter>" << mMovementParameter << "</MovementParameter> \n";

    // // Call direct parent class
    // AbstractForce<3>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////



// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(DriftPreventForce)
