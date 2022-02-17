#ifndef DriftPreventForce_HPP_
#define DriftPreventForce_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractOffLatticeCellPopulation.hpp"

/**
 * A force class to model random cell movement.
 */
class DriftPreventForce : public AbstractForce<3>
{
private :

    /**
     * Random Movement Parameter.
     */
    double mMovementParameter;

    double mTissueMiddle;

    /**
     * Archiving.
     */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<3> >(*this);
        archive & mMovementParameter;
        archive & mTissueMiddle;
    }

public :

    /**
     * Constructor.
     */
    DriftPreventForce();

    /**
     * Destructor.
     */
    ~DriftPreventForce();

    /**
     * Set the diffusion constant for the cells.
     *
     * @param movementParameter the movement parameter to use
     */
    void SetMovementParameter(double movementParameter);

    /**
     * Get the random motion coefficient.
     *
     * @return mMovementParameter
     */
    double GetMovementParameter();

    void SetTissueMiddle(double tissue_middle);

    double GetTissueMiddle();



    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the tissue
     *
     */
    void AddForceContribution(AbstractCellPopulation<3>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(DriftPreventForce)

#endif /*DriftPreventForce_HPP_*/
