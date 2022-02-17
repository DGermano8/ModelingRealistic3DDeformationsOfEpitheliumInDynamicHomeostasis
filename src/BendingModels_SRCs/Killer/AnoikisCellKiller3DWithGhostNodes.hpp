#ifndef ANOIKISCELLKILLER3DWITHGHOSTNODES_HPP_
#define ANOIKISCELLKILLER3DWITHGHOSTNODES_HPP_

#include "AbstractCellKiller.hpp"
#include "MeshBasedCellPopulation.hpp"
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
 * A cell killer that randomly kills cells based on the user set probability.
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
class AnoikisCellKiller3DWithGhostNodes : public AbstractCellKiller<3>
{
private:

    unsigned mCellsRemovedByAnoikis;
    double mCutOffLength;
    double mCellPopulationWidth;
    double mCellPopulationDepth;

    std::vector<c_vector<double,5> > mLocationsOfAnoikisCells;

    // The output file directory for the simulation data that corresponds to the number of cells
    // killed by anoikis
    out_stream mAnoikisOutputFile;

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

        archive & mCellsRemovedByAnoikis;
        archive & mCutOffLength;
        archive & mCellPopulationWidth;
        archive & mCellPopulationDepth;
//        archive & mXLocationsOfAnoikisCells;
        archive & mOutputDirectory;
    }

public:

    /**
     * Default constructor.
     *
     * @param pCellPopulation pointer to a tissue
     * @param sloughOrifice whether to slough compressed cells at crypt orifice
     */
    AnoikisCellKiller3DWithGhostNodes(AbstractCellPopulation<3>* pCellPopulation, double cut_off = 100.0, double tissueWidth = 10.0, double tissueDepth = 10.0);

	// Destructor
	~AnoikisCellKiller3DWithGhostNodes();

    void SetOutputDirectory(std::string outputDirectory);

    std::string GetOutputDirectory();

    std::set<unsigned> GetNeighbouringNodeIndices(unsigned nodeIndex);

    bool HasCellPoppedUp(unsigned nodeIndex);

    std::vector<c_vector<unsigned,2> > RemoveByAnoikis();

    /**
     *  Loops over and kills cells by anoikis
     */
    void CheckAndLabelCellsForApoptosisOrDeath();

    /* After each event of cell killing in CheckAndLabelCellsForApoptosisOrDeath(), the information of whether to kill each cell
     * or not is passed to this method which then increments the member variables corresponding to the total number of cells
     * killed by anoikis
     */
    void SetNumberCellsRemovedByAnoikis(std::vector<c_vector<unsigned,2> > cellsRemoved);

    /* Returns the total number of cells removed by anoikis
     */
    unsigned GetNumberCellsRemovedByAnoikis();

    /* Storing the locations of those epithelial cells that get removed by anoikis
     */
    void SetLocationsOfCellsRemovedByAnoikis(std::vector<c_vector<unsigned,2> > cellsRemoved);

    /* Returns the locations of those cells removed by anoikis (time -- node_index -- x -- y -- z)
     *
     */
    std::vector<c_vector<double,5> > GetLocationsOfCellsRemovedByAnoikis();

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
CHASTE_CLASS_EXPORT(AnoikisCellKiller3DWithGhostNodes)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a AnoikisCellKiller3DWithGhostNodes.
 */
template<class Archive >
inline void save_construct_data(
    Archive & ar, const AnoikisCellKiller3DWithGhostNodes * t, const unsigned int file_version)
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
    Archive & ar, AnoikisCellKiller3DWithGhostNodes * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<3>* p_tissue;
    ar >> p_tissue;

    // Invoke inplace constructor to initialise instance
    ::new(t)AnoikisCellKiller3DWithGhostNodes(p_tissue);
}
}
}

#endif /*ANOIKISCELLKILLER3DWITHGHOSTNODES_HPP_*/
