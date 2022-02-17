#ifndef PERIODICBOXBOUNDARYCONDITION3DGHOSTS_HPP_
#define PERIODICBOXBOUNDARYCONDITION3DGHOSTS_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A boundary condition class for use with 3d box simulations. This restricts cells
 * to the periodic domain, and pins stromal cells below a certain height.
 */

class PeriodicBoxBoundaryCondition3dGhosts : public AbstractCellPopulationBoundaryCondition<3>
{
private:

	/* Width of the periodic domain */
	double mCellPopulationWidth;
	
	/* Depth of the periodic domain */
	double mCellPopulationDepth;
	
	/* The maximum height below which to pin stromal cells */
	double mMaxHeightForPinnedCells;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<3> >(*this);
        archive & mCellPopulationWidth;
        archive & mCellPopulationDepth;
        archive & mMaxHeightForPinnedCells;
    }

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     */
	PeriodicBoxBoundaryCondition3dGhosts(AbstractCellPopulation<3>* pCellPopulation);

    /**
     * Destructor.
     */
    ~PeriodicBoxBoundaryCondition3dGhosts();

    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    void ImposeBoundaryCondition(const std::map<Node<3>*, c_vector<double, 3> >& rOldLocations);
    //void ImposeBoundaryCondition(const std::vector< c_vector<double, 3> >& rOldLocations);

    /**
     * Overridden VerifyBoundaryCondition() method.
     * Verify the boundary conditions have been applied.
     * This is called after ImposeBoundaryCondition() to ensure the condition is still satisfied.
     *
     * @return whether the boundary conditions are satisfied.
     */
    bool VerifyBoundaryCondition();

    /* Set method for the maximum height to which to pin stromal cells */
    void SetMaxHeightForPinnedCells(double maxHeightForPinnedCells);

    /* Get method for the maximum height to which to pin stromal cells */
    double GetMaxHeightForPinnedCells();
    
    /* Set method for the width of the periodic domain */
    void SetCellPopulationWidth(double cellPopulationWidth);
    
    /* Get method for the width of the periodic domain */
    double GetCellPopulationWidth();
    
    /* Set method for the depth of the periodic domain */
    void SetCellPopulationDepth(double cellPopulationDepth);
    
    /* Get method for the depth of the periodic domain */
    double GetCellPopulationDepth();

    /**
     * Overridden OutputCellPopulationBoundaryConditionParameters() method.
     * Output cell population boundary condition parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(PeriodicBoxBoundaryCondition3dGhosts)
//EXPORT_TEMPLATE_CLASS_SAME_DIMS(PeriodicBoxBoundaryCondition3d)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a PeriodicBoxBoundaryCondition3d.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const PeriodicBoxBoundaryCondition3dGhosts * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<3>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;
}

/**
 * De-serialize constructor parameters and initialize a PeriodicBoxBoundaryCondition3d.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, PeriodicBoxBoundaryCondition3dGhosts * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<3>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)PeriodicBoxBoundaryCondition3dGhosts(p_cell_population);
}
}
} // namespace ...

#endif /*PERIODICBOXBOUNDARYCONDITION3DGHOSTS_HPP_*/
