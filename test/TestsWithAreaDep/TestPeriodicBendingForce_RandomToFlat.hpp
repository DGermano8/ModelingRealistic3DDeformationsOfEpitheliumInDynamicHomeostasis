#ifndef TEST3DBOXMODEL_HPP_
#define TEST3DBOXMODEL_HPP_

#include <cxxtest/TestSuite.h>

#include "CellBasedSimulationArchiver.hpp"
#include "Timer.hpp"

#include "TrianglesMeshReader.hpp"
#include "OffLatticeSimulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "DomMeshBasedCellPopulationWithGhostNodes.hpp"
// #include "MeshBasedCellPopulationWithGhostNodes.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StromalCellMutationState.hpp"
#include "CellsGenerator.hpp"
#include "RandomMotionForce.hpp"
#include "PeriodicCryptModelInteractionForceWithGhostNodes.hpp"
// #include "PeriodicBendingForce3dHeightWithGhostNodes.hpp"
#include "PeriodicBendingForce3dHeightAreaCurvatureWithGhostNodes.hpp"
#include "SloughingCellKiller3DWithGhostNodes.hpp"
#include "AnoikisCellKiller3DWithGhostNodes.hpp"
#include "UniformCellKiller3dWithGhostNodes.hpp"
// #include "PeriodicBoxBoundaryCondition3d.hpp"
#include "PeriodicBoxBoundaryCondition3dGhosts.hpp"
#include "PeriodicStromalBoxBoundaryCondition3d.hpp"
#include "TrianglesMeshWriter.hpp"
#include "Debug.hpp"

#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "CellIdWriter.hpp"
#include "CellAngleWriter.hpp"
#include "CellAreaWriter.hpp"


#include "CellMutationStatesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellAncestor.hpp"

#include "CellPopulationEpithelialWriter.hpp"
#include "VoronoiDataWriter.hpp"

#include "PetscSetupAndFinalize.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "PlanarDivisionRule.hpp"
#include "DriftPreventForce.hpp"
#include "FixedEpithelialBoundary3d.hpp"

// #include "MeshModifier.hpp"
#include "MeshRemeshModifier.hpp"
#include "PeriodicRemeshCellsModifier.hpp"
#include "MeshRemeshCutoffModifier.hpp"
// #include "DensityDependantCellKillerStrip3DWithGhostNodesV2.hpp"
#include "DensityDependantCellKillerStrip3DWithGhostNodesSimpleV2.hpp"

#include "DomSimpleWntCellCycleModel.hpp"
#include "NodeVelocityWriter.hpp"

// Tests have been copied directly from TestOffLatticeSimulation3d.hpp - most likely that bits will get changed
// and broken by me along the way

class Test3dBoxModel : public AbstractCellBasedTestSuite
{
private:

public:
    
    /* A cube that consists of a block of stromal cells, and a single layer of epithelial cells.
     * Target curvature for the basement membrane is zero, and multiple division events occur.
     */
    void TestPeriodicCubeWithGhosts() throw (Exception)
    {
    	RandomNumberGenerator::Instance()->Reseed(10);

        std::vector<Node<3>*> nodes;

        std::string output_directory = "PeriodicBend_Random2Flat_10";

        unsigned width = 8;	   // x
        unsigned height = 10;      // y
        unsigned ghosts_bottom = 0;       // ghosts > depth
        unsigned ghosts_top = 1;       // ghosts > depth
        unsigned num_tissue_depth = 1;
        unsigned depth = num_tissue_depth + (ghosts_bottom + ghosts_top) + 1;        // z

        // Initialise the tissue in an equilibrum state
        double width_space = 01.0;
        double height_space = 01.0*sqrt(0.75);
        double ghost_sep = 1.0;
        double depth_space = 0.738431690356779*1.0; //Magic number for z-spaceing...
        unsigned cells_per_layer = width*height;
        unsigned cell_iter = 0;

        double periodic_width = (double) (width+0.0)*width_space;
        double periodic_height = (double) (height+0.0)*height_space;

        double  wnt_strip_width = 1.25;

        double tissue_base = 5.0; //Hieght of tissue to prevent drift
        double tissue_middle = 0.0;

        bool isGhost = false;
        double x_coordinate, y_coordinate, z_coordinate;

        double centre_x = 30;
        double centre_y = 30;

        double spring_strength = 20.0;
        double spring_cuttoff = 1.5;

        double radius =  0;//periodic_width+1.0;
        double target_curvature = -0.2; //maximum curvature is 0.2066 -> higher curvature means smaller sphere
        double beta_parameter = 0.2*spring_strength;
        double alpha_parameter = 1.01;

        double time_step = 0.001;
        double end_time = 2.0;
        double plot_step = 1.0;

        bool include_springs = true;
        bool include_bending = true;

        int is_transit[depth*height*width];
        int num_real_nodes = 0;

        double eps_z = 0.25;
        double eps_xy = 0.1;

        std::vector<unsigned> ghost_node_indices, real_node_indices;
        for (unsigned k=0; k<(depth + 1 ); k++) //+1 puts a layer of ghosts on the bottom
        {
            isGhost = false;
            if(k < ghosts_bottom || k > ghosts_bottom + num_tissue_depth)
            {
                isGhost = true;
            }
            if(k == depth)
            {
                isGhost = true;
            }


            for (unsigned j=0; j<height; j++)
            {    
                for (unsigned i=0; i<width; i++)
                {
                    is_transit[cell_iter] = 0;
                        
                    c_vector<double, 3> node_i_new_location;

                    x_coordinate = (double) (i + 0.5*(j%2 + k%2))*width_space    + eps_xy*(2*RandomNumberGenerator::Instance()->ranf() - 1);
                    y_coordinate = (double) j*height_space                       + eps_xy*(2*RandomNumberGenerator::Instance()->ranf() - 1);
                    
                    if( k == depth)
                    {
                        z_coordinate = (double) tissue_base + (-1.0)*depth_space + eps_z*(2*RandomNumberGenerator::Instance()->ranf() - 1);

                    }
                    else
                    {
                        z_coordinate = (double) tissue_base + k*depth_space      + eps_z*(2*RandomNumberGenerator::Instance()->ranf() - 1);
                    }    
                    if( pow(x_coordinate - 0.5*periodic_width,2)+ pow(y_coordinate - 0.5*periodic_height ,2) <= pow(1.0,2) )
                    {
                        is_transit[cell_iter] = 1;
                    }
                        
                    if(isGhost)
                    {
                        ghost_node_indices.push_back(cell_iter);
                    }
                    else
                    {
                        num_real_nodes++;
                        tissue_middle = tissue_middle + z_coordinate;
                        real_node_indices.push_back(cell_iter);
                    }
                    nodes.push_back(new Node<3>(cell_iter,  false,  x_coordinate, y_coordinate, z_coordinate));
                    cell_iter++;



                    }
                }
        }
        tissue_middle = tissue_middle/num_real_nodes;

        MutableMesh<3,3> mesh(nodes);
        //mesh.Translate(0.5, 0.5);

        std::vector<CellPtr> cells;

    	// First we sort the real node indices into increasing order (the last ones will correspond to the
    	// epithelial nodes)
    	sort(real_node_indices.begin(), real_node_indices.end());

        // Set up cells by iterating through the mesh nodes
        unsigned num_nodes = mesh.GetNumAllNodes();
        unsigned num_epithelial_cells = (width)*(height);
        unsigned num_tissue_cells = (width)*(height)*(num_tissue_depth);
        unsigned num_ghosts_bottom = (width)*(height)*(ghosts_bottom);
        unsigned num_ghosts_top = (width)*(height)*(ghosts_top);

        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_stromal_state(new StromalCellMutationState);

        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        // Initialise Tissue cells (Stromal)
		for (unsigned i=0; i<real_node_indices.size()-num_epithelial_cells; i++)
		{
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();

            CellPtr p_differentiated_cell(new Cell(p_stromal_state, p_model));
            p_differentiated_cell->SetCellProliferativeType(p_differentiated_type);
            p_model->SetDimension(3);

            double birth_time = - 10.0;
            p_differentiated_cell->SetBirthTime(birth_time);

			cells.push_back(p_differentiated_cell);
        }

        // Initialise Epithelial cells
        for (unsigned i=real_node_indices.size()-num_epithelial_cells; i<real_node_indices.size(); i++)
        {

        	FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            p_model->SetMaxTransitGenerations(1);
            p_model->SetSDuration(2); 
            p_model->SetG2Duration(2); 
            p_model->SetMDuration(2); 
            p_model->SetDimension(3);
            
            CellPtr p_epithelial_cell(new Cell(p_state, p_model));
            double birth_time = -8;

            // p_epithelial_cell->SetCellProliferativeType(p_transit_type);
            p_epithelial_cell->SetCellProliferativeType(p_differentiated_type);
			
            p_epithelial_cell->SetBirthTime(birth_time);
            p_epithelial_cell->SetApoptosisTime(1);
            
            p_epithelial_cell->InitialiseCellCycleModel();
            
            MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (i));
            p_epithelial_cell->SetAncestor(p_cell_ancestor);

            cells.push_back(p_epithelial_cell);

        }
        
        std::cout<< "number of cells comp = " << real_node_indices.size() << "\n";
        std::cout<< "number of cells all  = " << num_epithelial_cells+num_tissue_cells << "\n";
        std::cout<< "number of ghosts     = " << ghost_node_indices.size() << "\n";

        DomMeshBasedCellPopulationWithGhostNodes<3> cell_population(mesh, cells, real_node_indices); //ghost_sep
        cell_population.SetDampingConstantNormal(2.0);
        cell_population.SetDampingConstantMutant(2.0);


        // Set the division rule for our population to be the random direction division rule
        boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_division_rule_to_set(new PlanarDivisionRule());
        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

        OffLatticeSimulation<3> simulator(cell_population);

        // Pass an adaptive numerical method to the simulation
        boost::shared_ptr<AbstractNumericalMethod<3,3> > p_method(new ForwardEulerNumericalMethod<3,3>());
        p_method->SetUseAdaptiveTimestep(false);
        simulator.SetNumericalMethod(p_method);

        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAngleWriter>();
        cell_population.AddCellWriter<CellAreaWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();
        cell_population.AddPopulationWriter<CellPopulationEpithelialWriter>();
        cell_population.AddPopulationWriter<NodeVelocityWriter>();

        //To fix paraview
        cell_population.SetWriteVtkAsPointsDom(true);

        std::map<Node<3>*, c_vector<double,3> > node_locations_before;
        for (std::list<CellPtr>::iterator cell_iter = cell_population.rGetCells().begin();
             cell_iter != cell_population.rGetCells().end();
             ++cell_iter)
        {

            Node<3>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
            node_locations_before[p_node] = p_node->rGetLocation();
        }

        TRACE("Initialise cells");
        for (typename AbstractCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
         cell_iter != cell_population.End();
         ++cell_iter)
        {
            cell_iter->GetCellData()->SetItem("angle_curvature", 0.0);
            cell_iter->GetCellData()->SetItem("cell_area", 2.0*sqrt(3.0)/4.0);

        }
        TRACE("Done");

        MAKE_PTR_ARGS(PeriodicBoxBoundaryCondition3dGhosts, boundary_condition, (&cell_population));
        boundary_condition->SetCellPopulationWidth(periodic_width);
        boundary_condition->SetCellPopulationDepth(periodic_height);
        boundary_condition->SetMaxHeightForPinnedCells(0.0);       			      
        boundary_condition->ImposeBoundaryCondition(node_locations_before);
        simulator.AddCellPopulationBoundaryCondition(boundary_condition);

		// Create periodic spring force law
        MAKE_PTR(PeriodicCryptModelInteractionForceWithGhostNodes<3>, periodic_spring_force);
        periodic_spring_force->SetUseOneWaySprings(false); //turning this on makes the stromal cells act as ghosts..
        periodic_spring_force->SetCutOffLength(spring_cuttoff);
        //                     SetEpithelialStromalCellDependentSprings(ind , Ep-Ep, Str-Str, Ep-Str, apcTwoHitStromalMultiplier);
        periodic_spring_force->SetEpithelialStromalCellDependentSprings(true, 1.0,     1.0,     0.5,    1.0);
        periodic_spring_force->SetPeriodicDomainWidth(periodic_width);
        periodic_spring_force->SetPeriodicDomainDepth(periodic_height);
        periodic_spring_force->SetMeinekeSpringStiffness(spring_strength);
        if(include_springs)
        {
            simulator.AddForce(periodic_spring_force);
        }

		// Create periodic basement membrane force law
        MAKE_PTR(PeriodicBendingForce3dHeightAreaCurvatureWithGhostNodes, periodic_bending_force);
        periodic_bending_force->SetOutputDirectory(output_directory);
        periodic_bending_force->SetHeightDependantCurvatureParameter(1.2);
        periodic_bending_force->SetBasementMembraneParameter(beta_parameter);
        periodic_bending_force->SetExponentParameter(alpha_parameter);
        // periodic_bending_force->SetCircularNonZeroTargetCurvatureRegion(true, target_curvature, radius, 20, 20);
        periodic_bending_force->SetCircularNonZeroTargetCurvatureRegion(true, target_curvature, radius, centre_x, centre_y);
        periodic_bending_force->SetPeriodicDomainWidth(periodic_width);
        periodic_bending_force->SetPeriodicDomainDepth(periodic_height);
        if(include_bending)
        {
            simulator.AddForce(periodic_bending_force);
        }

        simulator.SetOutputDirectory(output_directory);	 

        MAKE_PTR(MeshRemeshCutoffModifier<3>, pc_modifier);
        pc_modifier->SetOutputDirectory(output_directory + "/results_from_time_0");
        pc_modifier->SetWidth(periodic_width);
        pc_modifier->SetDepth(periodic_height);
        pc_modifier->SetCutoff(spring_cuttoff+0.001);
        simulator.AddSimulationModifier(pc_modifier);

        TRACE("Set uo done - lets solve");
    	
        simulator.SetSamplingTimestepMultiple(plot_step);			// Every hour
		simulator.SetEndTime(end_time);
        simulator.SetDt(time_step);

        Timer::Reset();
        simulator.Solve();
        Timer::Print("Time Ellapsed");
        
    }     

};

#endif /*TEST3DBOXMODEL_HPP_*/

