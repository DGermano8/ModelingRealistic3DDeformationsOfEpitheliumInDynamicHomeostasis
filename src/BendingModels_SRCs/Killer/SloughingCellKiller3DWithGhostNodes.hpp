#ifndef SLOUGHINGCELLKILLER3DWITHGHOSTNODES_HPP_
#define SLOUGHINGCELLKILLER3DWITHGHOSTNODES_HPP_

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
 * A cell killer that randomly kills cells that move beyond the edges of the tissue box
 */
class SloughingCellKiller3DWithGhostNodes : public AbstractCellKiller<3>
{
private:

	/* x and y boundaries for sloughing epithelial cells that move beyond the edges of the tissue box */
    double mCellPopulationWidth;
	double mCellPopulationDepth;

    unsigned mCellsRemovedBySloughing;

    std::vector<c_vector<double,5> > mLocationsOfSloughedCells;

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

        archive & mCellPopulationWidth;
        archive & mCellPopulationDepth;
        archive & mCellsRemovedBySloughing;
//        archive & mLocationsOfSloughedCells;
        archive & mOutputDirectory;
    }

public:

    /**
     * Default constructor.
     *
     * @param pCellPopulation pointer to a tissue
     * @param sloughOrifice whether to slough compressed cells at crypt orifice
     */
    SloughingCellKiller3DWithGhostNodes(AbstractCellPopulation<3>* pCellPopulation, double tissueWidth = 10.0, double tissueDepth = 10.0);

	// Destructor
	~SloughingCellKiller3DWithGhostNodes();

    double GetCellPopulationWidth() const;

    double GetCellPopulationDepth() const;

    void SetOutputDirectory(std::string outputDirectory);

    std::string GetOutputDirectory();

    std::vector<c_vector<unsigned,2> > RemoveBySloughing();

    /**
     *  Loops over and kills cells by sloughing
     */
    void CheckAndLabelCellsForApoptosisOrDeath();

    /* After each event of cell killing in CheckAndLabelCellsForApoptosisOrDeath(), the information of whether to kill each cell
     * or not is passed to this method which then increments the member variables corresponding to the total number of cells
     * killed by sloughing
     */
    void SetNumberCellsRemovedBySloughing(std::vector<c_vector<unsigned,2> > cellsRemoved);

    /* Returns the total number of cells removed by sloughing
     *
     */
    unsigned GetNumberCellsRemovedBySloughing();

    /* Storing the locations of those epithelial cells that get removed by sloughing
     *
     */
    void SetLocationsOfCellsRemovedBySloughing(std::vector<c_vector<unsigned,2> > cellsRemoved);

    /* Returns the x-coordinates of those cells removed by anoikis
     *
     */
    std::vector<c_vector<double,5> > GetLocationsOfCellsRemovedBySloughing();

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
CHASTE_CLASS_EXPORT(SloughingCellKiller3DWithGhostNodes)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a SloughingCellKiller3DWithGhostNodes.
 */
template<class Archive >
inline void save_construct_data(
    Archive & ar, const SloughingCellKiller3DWithGhostNodes * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<3>* const p_tissue = t->GetCellPopulation();
    ar << p_tissue;
    double tissue_width = t->GetCellPopulationWidth();
    ar << tissue_width;
    double tissue_depth = t->GetCellPopulationDepth();
    ar << tissue_depth;
}

/**
 * De-serialize constructor parameters and initialise Crypt.
 */
template<class Archive >
inline void load_construct_data(
    Archive & ar, SloughingCellKiller3DWithGhostNodes * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<3>* p_tissue;
    ar >> p_tissue;
    double tissue_width;
    ar >> tissue_width;
    double tissue_depth;
    ar >> tissue_depth;

    // Invoke inplace constructor to initialise instance
    ::new(t)SloughingCellKiller3DWithGhostNodes(p_tissue, tissue_width, tissue_depth);
}
}
}

#endif /*SLOUGHINGCELLKILLER3DWITHGHOSTNODES_HPP_*/
