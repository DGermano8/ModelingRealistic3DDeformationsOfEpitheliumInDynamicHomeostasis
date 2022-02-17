#ifndef UNIFORMCELLKILLER3DWITHGHOSTNODES_HPP_
#define UNIFORMCELLKILLER3DWITHGHOSTNODES_HPP_

#include "AbstractCellKiller.hpp"
#include "DomMeshBasedCellPopulationWithGhostNodes.hpp"
#include "OutputFileHandler.hpp"
#include "StromalCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApoptoticCellProperty.hpp"
#include "SimpleWntCellCycleModel.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A cell killer that randomly kills cells in the outer rim of
 * the epithelial monolayer based on the user set probability.
 *
 * The probability passed into the constructor will be the probability
 * of any cell dying whenever CheckAndLabelCellsForApoptosis() is called.
 *
 * Note this does take into account timesteps - the input probability is the
 * probability that in an hour's worth of trying, the cell killer will have
 * successfully killed a given cell. In the method CheckAndLabelSingleCellForApoptosis()
 * this probability is used to calculate the probability that the cell is killed
 * at a given time step.
 */
class UniformCellKiller3dWithGhostNodes : public AbstractCellKiller<3>
{
private:

    /**
     * Probability that in an hour's worth of trying, this cell killer
     * will have successfully killed a given cell.
     */
    double mProbabilityOfDeathInAnHour;

    /* Boundary in x-direction below which (i.e. with a lower x coord) cells are randomly killed */
    double mMinXBoundary;

    /* Boundary in x-direction beyond which (i.e. with a greater x coord) cells are randomly killed */
    double mMaxXBoundary;

    /* Boundary in x-direction below which (i.e. with a lower x coord) cells are randomly killed */
    double mMinYBoundary;

    /* Boundary in y-direction beyond which (i.e. with a greater y coord) cells are randomly killed */
    double mMaxYBoundary;

    /* minimum cell population to know when to stop killing */
    double mMinCellPopulation;

    // The output file directory for the simulation data that corresponds to the number of cells
    // killed by anoikis
    out_stream mRandomDeathOutputFile;

    std::string mOutputDirectory;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * Serialization of singleton objects must be done with care.
     * Before the object is serialized via a pointer, it *MUST* be
     * serialized directly, or an assertion will trip when a second
     * instance of the class is created on de-serialization.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<3> >(*this);

        archive & mProbabilityOfDeathInAnHour;
        archive & mMinXBoundary;
        archive & mMaxXBoundary;
        archive & mMinYBoundary;
        archive & mMaxYBoundary;
        archive & mMinCellPopulation;
        archive & mOutputDirectory;
    }

public:

    /**
     * Default constructor.
     *
     * @param pCellPopulation pointer to a tissue
     * @param sloughOrifice whether to slough compressed cells at crypt orifice
     */
    UniformCellKiller3dWithGhostNodes(AbstractCellPopulation<3>* pCellPopulation, double probabilityOfDeathInAnHour = 0.5,
    		double minXBoundary = 0.0, double maxXBoundary = 50.0, double minYBoundary = 0.0, double maxYBoundary = 50.0, double cellPop = 0.0);

	// Destructor
	~UniformCellKiller3dWithGhostNodes();


    /*@ return mProbabilityOfDeathInAnHour */
    double GetDeathProbabilityInAnHour() const;

    /* Get max and min boundaries in x and y directions */
    double GetXMinBoundary() const;

    double GetXMaxBoundary() const;

    double GetYMinBoundary() const;

    double GetYMaxBoundary() const;

    void SetOutputDirectory(std::string outputDirectory);

    std::string GetOutputDirectory();

    std::vector<c_vector<unsigned,2> > RemoveByRandomSelection();

    /**
     * Overridden method to test a given cell for apoptosis.
     *
     * @param pCell the cell to test for apoptosis
     */
    void CheckAndLabelSingleCellForApoptosis(CellPtr pCell);

    /**
     *  Loops over and kills cells by anoikis
     */
    void CheckAndLabelCellsForApoptosisOrDeath();

    /**
     * Outputs cell killer parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellKillerParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(UniformCellKiller3dWithGhostNodes)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a RandomCellKiller3d.
 */
template<class Archive >
inline void save_construct_data(
    Archive & ar, const UniformCellKiller3dWithGhostNodes * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<3>* const p_tissue = t->GetCellPopulation();
    ar << p_tissue;
}

/**
 * De-serialize constructor parameters and initialise Crypt.
 */
template<class Archive >
inline void load_construct_data(
    Archive & ar, UniformCellKiller3dWithGhostNodes * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<3>* p_tissue;
    ar >> p_tissue;

    // Invoke inplace constructor to initialise instance
    ::new(t)UniformCellKiller3dWithGhostNodes(p_tissue);
}
}
}

#endif /* UNIFORMCELLKILLER3DWITHGHOSTNODES_HPP_ */
