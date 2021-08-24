///////////////////////////////////////////////////////////////////////////////
//
// File: TimeIntegrationSchemeFIT.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Header file of time integration scheme FIT base class
//
///////////////////////////////////////////////////////////////////////////////

// Note: The file is named TimeIntegrationSchemeFIT to parallel the
// TimeIntegrationSchemeGLM file but the class is named
// FractionalInTimeIntegrationScheme so keep with the factory naming
// convention.

#ifndef NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_TIME_INTEGRATION_SCHEME_FIT
#define NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_TIME_INTEGRATION_SCHEME_FIT

#define LUE LIB_UTILITIES_EXPORT

#include <string>

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>

namespace Nektar
{
namespace LibUtilities
{

///////////////////////////////////////////////////////////////////////////////
/// Class for fractional-in-time integration.
class FractionalInTimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    /// Constructor
    FractionalInTimeIntegrationScheme(std::string variant, unsigned int order,
                                      std::vector<NekDouble> freeParams);

    /// Destructor
    virtual ~FractionalInTimeIntegrationScheme()
    {
    }

    /// Creator
    static TimeIntegrationSchemeSharedPtr
    create(std::string variant,
           unsigned int order,
           std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            FractionalInTimeIntegrationScheme>::AllocateSharedPtr(variant,
                                                                  order,
                                                                  freeParams);

        return p;
    }

    static std::string className;

    // Access methods from the base class that are virtual
    LUE virtual std::string GetName() const
    {
        return m_name;
    }

    LUE virtual std::string GetVariant() const
    {
        return m_variant;
    }

    LUE virtual unsigned int GetOrder() const
    {
        return m_order;
    }

    LUE virtual std::vector< NekDouble > GetFreeParams() const
    {
        return m_freeParams;
    }

    LUE virtual TimeIntegrationSchemeType GetIntegrationSchemeType() const
    {
        return m_schemeType;
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE unsigned int GetNumIntegrationPhases() const
    {
        return 1;
    }

    /**
     * \brief Gets the solution vector of the ODE
     */
    inline const TripleArray &GetSolutionVector() const
    {
        return m_u;
    }

    /**
     * \brief Sets the solution vector of the ODE
     */
    inline void SetSolutionVector(const int Offset, const DoubleArray &y)
    {
        m_u[Offset] = y;
    }

    // The worker methods from the base class that are virtual
    LUE virtual void InitializeScheme(
        const NekDouble deltaT, ConstDoubleArray &y_0,
        const NekDouble time, const TimeIntegrationSchemeOperators &op);

    LUE virtual ConstDoubleArray &TimeIntegrate(
        const int timestep, const NekDouble delta_t,
        const TimeIntegrationSchemeOperators &op);

    LUE virtual void print(std::ostream &os) const;
    LUE virtual void printFull(std::ostream &os) const;

    // Friend classes
    LUE friend std::ostream &operator<<(std::ostream &os,
        const FractionalInTimeIntegrationScheme &rhs);
    LUE friend std::ostream &operator<<(std::ostream &os,
        const FractionalInTimeIntegrationSchemeSharedPtr &rhs);

protected:
    struct Instance
    {
        int base;

        int  index;           // Index of this instance
        bool active;          // Used to determine if active
        int  activecounter;   // counter used to flip active bit
        int  activebase;

        // Major storage for auxilliary ODE solutions.
        // Storage for values of y currently used to update u
        ComplexTripleArray    stage_y;
        std::pair< int, int > stage_ind;  // Time-step counters indicating the
                                          // interval ymain is associated with

        // Staging allocation
        bool stage_active;
        int stage_ccounter;
        int stage_cbase;    // This base is halved after the first cycle
        int stage_fcounter;
        int stage_fbase;    // This base is halved after the first cycle

        // Ceiling stash allocation
        int                   cstash_counter;  // Counter used to determine
                                               // when to stash
        int                   cstash_base;     // base for counter
        ComplexTripleArray    cstash_y;
        std::pair< int, int > cstash_ind;      // ind(1) is never used:
                                               // it always matches main.ind(1)

        // Ceiling sandbox allocation
        bool                  csandbox_active; // Flag to determine when
                                               // stash 2 is utilized
        int                   csandbox_counter;
        ComplexTripleArray    csandbox_y;
        std::pair< int, int > csandbox_ind;

        // Floor stash
        int                   fstash_base;
        ComplexTripleArray    fstash_y;
        std::pair< int, int > fstash_ind;

        // Floor sandbox
        bool                  fsandbox_active;
        int                   fsandbox_activebase;
        int                   fsandbox_stashincrement;
        ComplexTripleArray    fsandbox_y;
        std::pair< int, int > fsandbox_ind;

        // Talbot quadrature rule
        ComplexSingleArray z;
        ComplexSingleArray w;

        TripleArray As;

        ComplexSingleArray E;
        ComplexDoubleArray Eh;
        ComplexDoubleArray AtEh;
    };

    inline unsigned int modIncrement(const unsigned int counter,
                                     const unsigned int base) const;

    inline unsigned int computeL( const unsigned int base,
                                  const unsigned int m ) const;

    inline unsigned int  computeQML( const unsigned int base,
                                     const unsigned int m );

    inline unsigned int computeTaus( const unsigned int base,
                                     const unsigned int m );

    void talbotQuadrature(const unsigned int nQuadPts,
                          const NekDouble mu,
                          const NekDouble nu,
                          const NekDouble sigma,
                                ComplexSingleArray &lamb,
                                ComplexSingleArray &w) const;

    void integralClassInitialize(const unsigned int index,
                                 Instance &instance) const;

    void updateStage(const unsigned int timeStep,
                           Instance &instance);

    void finalIncrement(const unsigned int timeStep,
                        const TimeIntegrationSchemeOperators &op);

    void integralContribution(const unsigned int timeStep,
                              const unsigned int tauml,
                              const Instance &instance);

    void timeAdvance(const unsigned int timeStep,
                     const TimeIntegrationSchemeOperators &op,
                           Instance &instance,
                           ComplexTripleArray &y);

    void advanceSandbox(const unsigned int timeStep,
                        const TimeIntegrationSchemeOperators &op,
                              Instance &instance);

    // Variables common to all schemes.
    std::string m_name;
    std::string m_variant;
    unsigned int m_order{0};
    std::vector< NekDouble > m_freeParams;

    TimeIntegrationSchemeType m_schemeType{eFractionalInTime};

    // Varaibles and methods specific to FIT integration schemes.
    NekDouble    m_deltaT{0};
    NekDouble    m_T{0};          // Finial time
    unsigned int m_maxTimeSteps;  // Number of time steps.
    NekDouble    m_alpha{0.3};    // Value for exp integration.
    unsigned int m_base{4};       // "Base" of the algorithm.
    unsigned int m_nQuadPts{20};  // Number of Talbot quadrature rule points
    NekDouble    m_sigma{0};
    NekDouble    m_mu0{8};
    NekDouble    m_nu{0.6};

    int m_nvars{0};            // Number of variables in the integration scheme.
    int m_npoints{0};          // Number of points    in the integration scheme.

    unsigned int m_Lmax{0};    // Maxium number of integral groups.
    Array<OneD, Instance> m_integral_classes;

    // Demarcation integers
    Array<OneD, int> m_qml;
    // Demarcation interval markers
    Array<OneD, int> m_taus;

    // Storage of the initial values.
    DoubleArray m_u0;
    // Storage of the next solution from the final increment.
    DoubleArray m_uNext;
    // Storage for the integral contribution.
    ComplexDoubleArray m_uInt;
    // Storage for the exponential factor in the integral contribution.
    ComplexSingleArray m_expFactor;

    // Storage of previous states and associated timesteps.
    TripleArray m_u;

    // Storage for the stage derivative as the data will be re-used to
    // update the solution.
    DoubleArray m_F;

    // J
    SingleArray m_J;

    // Ahat array one for each order.
    TripleArray m_Ahats;

    // Mulitply the last Ahat array, transposed by J
    SingleArray m_AhattJ;

}; // end class FractionalInTimeIntegrator

LUE std::ostream &operator<<(
    std::ostream &os,
    const FractionalInTimeIntegrationScheme &rhs);
LUE std::ostream &operator<<(
    std::ostream &os,
    const FractionalInTimeIntegrationSchemeSharedPtr &rhs);

} // end namespace LibUtilities
} // end namespace Nektar

#endif
