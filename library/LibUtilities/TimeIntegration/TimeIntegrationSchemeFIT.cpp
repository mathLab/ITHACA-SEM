///////////////////////////////////////////////////////////////////////////////
//
// File FractionalInTimeIntegrationScheme.cpp
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
// Description: implementation of time integration key class
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeFIT.h>

namespace Nektar
{
namespace LibUtilities
{

std::ostream &operator<<(std::ostream &os, const FractionalInTimeIntegrationScheme &rhs)
{
    os << "Time Integration Scheme: " << rhs.GetFullName() << std::endl
       << "Base " << rhs.m_base << std::endl
       << "Number of quadature points " << rhs.m_nQuadPts << std::endl
       << "Alpha " << rhs.m_alpha << std::endl;

    return os;
}

std::ostream &operator<<(std::ostream &os,
                         const FractionalInTimeIntegrationSchemeSharedPtr &rhs)
{
    os << *rhs.get();

    return os;
}

void FractionalInTimeIntegrationScheme::
InitializeScheme(const NekDouble deltaT,
                       ConstDoubleArray &y_0,
                 const NekDouble time,
                 const TimeIntegrationSchemeOperators &op)

{
    boost::ignore_unused(op);

    m_nvars   = y_0.num_elements();
    m_npoints = y_0[0].num_elements();

    m_deltaT = deltaT;

    m_T = time;  // Finial time;
    m_maxTimeSteps = m_T / m_deltaT;

    // The +2 below is a demonstration of how to do additional
    // allocation in case one decides later to increase T
    m_Lmax = compute_L( m_base, m_maxTimeSteps ) + 2;

    // Demarcation integers - one array that is re-used
    m_qml  = Array<OneD, int>( m_Lmax-1, 0 );
    // Demarcation interval markers - one array that is re-used
    m_taus = Array<OneD, int>( m_Lmax+1, 0 );

    // Storage of the initial values.
    m_u0    = Array<TwoD,              NekDouble>  (m_nvars, m_npoints, 0.0);

    for( int i=0; i<m_nvars; ++i )
    {
        for( int j=0; j<m_npoints; ++j )
        {
            m_u0[i][j] = y_0[i][j];
        }
    }

    // Storage of the next solution from the final increment.
    m_uNext = Array<TwoD,              NekDouble>  (m_nvars, m_npoints, 0.0);
    // Storage for the integral contribution.
    m_uInt  = Array<TwoD, std::complex<NekDouble> >(m_nvars, m_npoints, 0.0);

    // Storage for the exponential factor in the integral
    // contribution. One array that is re-used
    m_expFactor = Array<OneD, std::complex<NekDouble> >(m_nQuadPts, 0.0);

    // Storage of the derivatives. One array that is re-used
    m_F = TripleArray( m_order+1 );

    // Storage of previous states and associated timesteps.
    m_u = TripleArray( m_order+1 );

    for( unsigned int m=0; m<=m_order; ++m )
    {
        m_F[m] = DoubleArray( m_nvars );
        m_u[m] = DoubleArray( m_nvars );

        for( unsigned int i=0; i<m_nvars; ++i )
        {
            m_F[m][i] = SingleArray( m_npoints, 0.0 );
            m_u[m][i] = SingleArray( m_npoints, 0.0 );
        }
    }

    // J
    m_J = Array<OneD, NekDouble>(m_order, 0.0);

    m_J[0] = pow( m_deltaT, m_alpha ) / tgamma( m_alpha+1. );

    for( unsigned int m=1, m_1=0; m<m_order; ++m, ++m_1 )
    {
        m_J[m] = m_J[m_1] * NekDouble(m) / (NekDouble(m) + m_alpha);
    }

    // Ahat array, one for each order.
    m_Ahats = Array<OneD, Array<TwoD, NekDouble>>(m_order+1);

    Array<OneD, Array<TwoD, NekDouble>> &Ahats = m_Ahats;

    for( unsigned int m=1; m<=m_order; ++m )
    {
        Ahats[m] = Array<TwoD, NekDouble>(m, m, 0.0);

        switch( m )
        {
          case 1:
            Ahats[m][0][0] = 1.;
            break;

          case 2:
            Ahats[m][0][0] =  1.;      Ahats[m][0][1] =  0.;
            Ahats[m][1][0] =  1.;      Ahats[m][1][1] = -1.;
            break;

          case 3:
            Ahats[m][0][0] =  1.;     Ahats[m][0][1] =  0.;  Ahats[m][0][2] = 0;
            Ahats[m][1][0] =  3./2.;  Ahats[m][1][1] = -2.;  Ahats[m][1][2] = 1./2.;
            Ahats[m][2][0] =  1./2.;  Ahats[m][2][1] = -1.;  Ahats[m][2][2] = 1./2.;
            break;

          case 4:
            Ahats[m][0][0] =  1.;      Ahats[m][0][1] =  0.;
            Ahats[m][0][2] =  0.;      Ahats[m][0][3] =  0.;

            Ahats[m][1][0] =  11./6.;  Ahats[m][1][1] = -3;
            Ahats[m][1][2] =  3./2.;   Ahats[m][1][3] = -1./3.;

            Ahats[m][2][0] =  1.;      Ahats[m][2][1] = -5./2.;
            Ahats[m][2][2] =  2.;      Ahats[m][2][3] = -1./2.;

            Ahats[m][3][0] =  1./6.;   Ahats[m][3][1] = -1./2.;
            Ahats[m][3][2] =  1./2.;   Ahats[m][3][3] = -1./6.;
            break;

          default:

            Ahats[m][0][0] = 1;

            for( unsigned int j=2; j<=m; ++j )
            {
                for( unsigned int i=0; i<m; ++i )
                {
                    Ahats[m][j-1][i] = pow( (1-j), i );
                }
            }

            ASSERTL1(false, "No matrix inverse.");

            // Ahats[m] = inv(Ahats[m]);

            break;
        }
    }

    // Mulitply the last Ahat array, transposed, by J
    m_AhattJ = Array<OneD, NekDouble>(m_order, 0.0);

    for( unsigned int i=0; i<m_order; ++i )
    {
        for( unsigned int j=0; j<m_order; ++j )
        {
            m_AhattJ[i] += Ahats[m_order][j][i] * m_J[j];
        }
    }

    m_integral_classes = Array<OneD, Instance>(m_Lmax);

    for (int l=0; l<m_Lmax; ++l)
    {
        integral_class_initialize( l+1, m_integral_classes[l] );
        // std::cout << "initialize_integral_class out integral_classes[i] "
        //            << m_integral_classes[l].index << std::endl;

        // for( int i=0; i<m_nvars; ++i )
        //   {
        //      for( int j=0; j<m_nQuadPts/4; ++j )
        //        {
        //          std::cout << m_integral_classes[l].stage_y[i][j] << "  ";
        //        }
        //      std::cout << std::endl;
        //   }
        // std::cout << std::endl;
    }
}


ConstDoubleArray &
FractionalInTimeIntegrationScheme::
TimeIntegrate(const int timestep, const NekDouble delta_t,
              const TimeIntegrationSchemeOperators &op)
{
    boost::ignore_unused(delta_t);

    ASSERTL1( delta_t == m_deltaT,
              "Delta T has changed which is not permitted." );

    // The Fractional in Time works via the logical? time step value.
    int timeStep = timestep + 1;

    // std::cout << "order  " << m_order << "  "
    //        << "timeStep  " << timeStep << std::endl;

    // Update the storage and counters for integral classes.  Performs
    // staging for updating u.
    for (int l=0; l<m_Lmax; ++l)
    {
        // std::cout << "update_stage in m_integral_classes[l] "
        //          << m_integral_classes[l].index << std::endl;

        // for( int i=0; i<m_nvars; ++i )
        //   {
        //    for( int j=0; j<m_nQuadPts/4; ++j )
        //      {
        //        std::cout << m_integral_classes[l].stage_y[i][j] << "  ";
        //      }
        //    std::cout << std::endl;
        //   }
        // std::cout << std::endl;

        update_stage(timeStep, m_integral_classes[l]);

        // std::cout << "update_stage out m_integral_classes[i] "
        //          << m_integral_classes[l].index << std::endl;
        // for( int i=0; i<m_nvars; ++i )
        //   {
        //    for( int j=0; j<m_nQuadPts/4; ++j )
        //      {
        //        std::cout << m_integral_classes[l].stage_y[i][j] << "  ";
        //      }
        //    std::cout << std::endl;
        //   }
        // std::cout << std::endl;
    }

    // Compute u update to time timeStep * m_deltaT.  Stored in
    // m_uNext.
    final_increment(timeStep, op);

    // for( int i=0; i<m_nvars; ++i )
    //   {
    //  for( int j=0; j<m_npoints; ++j )
    //    {
    //      std::cout << m_uNext[i] << "  ";
    //    }
    //  std::cout << std::endl;
    //   }
    // std::cout << std::endl;

    // Contributions to the current integral
    int L = compute_taus( m_base, timeStep );

    for (int l=0; l<L; ++l)
    {
        // Integral contribution over [taus(i+1) taus(i)]. Stored in
        // m_uInt.
        integral_contribution( timeStep, m_taus[l], m_integral_classes[l] );

        // for( int i=0; i<m_nvars; ++i )
        //   {
        //    for( int j=0; j<m_npoints; ++j )
        //      {
        //        std::cout << m_uInt[i][j] << "  ";
        //      }
        //    std::cout << std::endl;
        //   }
        // std::cout << std::endl;

        for( int i=0; i<m_nvars; ++i )
        {
            for( int j=0; j<m_npoints; ++j )
            {
                m_uNext[i][j] += m_uInt[i][j].real();
            }
        }
    }

    // Shuffle the previous solutions back one in the history.
    for( int m=m_order; m>0; --m )
    {
        for( int i=0; i<m_nvars; ++i )
        {
            for( int j=0; j<m_npoints; ++j )
            {
                m_u[m][i][j] = m_u[m-1][i][j];
            }
        }
    }

    // Get the current solution.
    for( int i=0; i<m_nvars; ++i )
    {
        for( int j=0; j<m_npoints; ++j )
        {
            m_u[0][i][j] = m_uNext[i][j] + m_u0[i][j];
        }
    }

    // Dump the current solution.
    // for( int j=0; j<m_npoints; ++j )
    // {
    //     for( int i=0; i<m_nvars; ++i )
    //     {
    //         for( int m=0; m<m_order+1; ++m )
    //         {
    //             std::cout << m_u[m][i][j] << "  ";
    //         }
    //         std::cout << std::endl;
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    // Update the storage and counters for integral classes to
    // time timeStep * m_deltaT. Also update the sandboxes and
    // stashes.
    for (int i=0; i<m_Lmax; ++i)
    {
        advance_sandbox( timeStep, op, m_integral_classes[i] );
    }

    return m_u[0];
}


unsigned int FractionalInTimeIntegrationScheme::
modincrement(const int unsigned counter,
             const int unsigned base) const
{
   return (counter+1) % base;
}

// Computes the smallest integer L such that base < 2 * base^l.
unsigned int FractionalInTimeIntegrationScheme::
compute_L( const unsigned int base,
           const unsigned int l ) const
{
    unsigned int L = ceil(log(l/2.0) / log(base));

    if( l % (unsigned int)(2 * pow(base,L)) == 0 )
    {
        ++L;
    }

    return L;
}

//  Computes the demarcation integers q_{m, ell}.
//
//  Returns a length-(L-1) vector qml such that h*taus are interval
//  boundaries for a partition of [0, m h]. The value of h is not
//  needed to compute this vector.
unsigned int
FractionalInTimeIntegrationScheme::compute_qml( const unsigned int base,
                                                    const unsigned int m )
{
    int L = compute_L(base, m);

    // qml = Array<OneD, int>( L-1, 0 );
    
    for( unsigned int i=0; i<L-1; ++i )
    {
        m_qml[i] = floor(m / pow(base, i+1) ) - 1;

        // std::cout << qml[i] << "  ";
    }
    // std::cout << std::endl;
    
    return L;
}

// Computes the demarcation interval marker tau_{m, ell}.
//
// Returns a length-(L+1) vector tauml such that h*taus are interval
// boundaries for a partition of [0, m h]. The value of h is not
// needed to compute this vector.
unsigned int
FractionalInTimeIntegrationScheme::compute_taus( const unsigned int base,
                                                     const unsigned int m )
{
    if( m == 1 )
    {
        m_taus[0] = 0;

        return 0;
    }
    else
    {
        unsigned int L = compute_qml(base, m);

        // taus = Array<OneD, int>( L+1, 0 );

        m_taus[0] = m - 1;

        for( unsigned int i=1; i<L; ++i )
        {
            m_taus[i] = m_qml[i-1] * pow(base, i);
        }

        m_taus[L] = 0;

        // std::cout << "taus " << L << std::endl;
        // for( unsigned int i=0; i<=L; ++i )
        //   std::cout << taus[i] << std::endl;
        // std::cout << std::endl;

        return L;
    }
}


// talbot_quadrature -- returns quadrature rule over Tablot contour
//
// Returns a quadrature rule over the Tablot contour defined by the
// parameterization.
//
// gamma(th) = sigma + mu * ( th*cot(th) + i*nu*th ),  -pi < th < pi
//
// An N-point rule is returned, equidistant in the parameter
// theta. Denoting the contour Gamma, the returned quadrature rule is
// the approximation
//
void FractionalInTimeIntegrationScheme::
talbot_quadrature(const unsigned int nQuadPts,
                  const NekDouble mu,
                  const NekDouble nu,
                  const NekDouble sigma,
                  Array<OneD, std::complex<NekDouble> > &lamb,
                  Array<OneD, std::complex<NekDouble> > &w) const
{
    lamb = Array<OneD, std::complex<NekDouble> >(nQuadPts, 0.0);
    w    = Array<OneD, std::complex<NekDouble> >(nQuadPts, 0.0);

    for( unsigned int q=0; q<nQuadPts; ++q )
    {
        NekDouble th =
          (NekDouble(q) + 0.5) / NekDouble(nQuadPts) * 2.0 * M_PI - M_PI;

        lamb[q] = sigma + mu * th * std::complex<NekDouble>(1./tan(th), nu);

        w[q] = std::complex<NekDouble>(0, -1./NekDouble(nQuadPts)) *
          mu * std::complex<NekDouble>(1./tan(th) - th/(sin(th)*sin(th)), nu);
    }

    // Special case for th = 0 which happens when there is an odd number
    // of quadrature points.
    if( nQuadPts % 2 == 1 )
    {
        unsigned int q = (nQuadPts+1) / 2;

        lamb[q] = std::complex<NekDouble>(sigma + mu, 0);

        w[q] = std::complex<NekDouble>(nu * mu / nQuadPts, 0);
    }
}


void FractionalInTimeIntegrationScheme::
integral_class_initialize(const unsigned int index,
                          Instance &instance) const
{
    // This object stores information for performing integration over an
    // interval [a, b]. (Defined by taus in the parent calling
    // function.)

    // The "main" object stores information about [a,b]. In particular,
    // main.ind identifies [a,b] via multiples of h.

    // Periodically the values of [a,b] need to be incremented. The
    // necessary background storage to accomplish this increment depends
    // whether a or b is being incremented.

    // The objects with "f" ("Floor") modifiers are associated with
    // increments of the interval floor a.

    // The objects with "c" ("Ceiling") modifiers are associated with
    // increments of the interval ceiling b.

    // Items on the "stage" are stored for use in computing u at the
    // current time.  Items in the "stash" are stored for use for future
    // staging Items in the "sandbox" are being actively updated at the
    // current time for future stashing.

    // This is the same for all integral classes, so there's probably a
    // better way to engineer this. And technically, all that's needed
    // is the array K(instance.z) anyway.
    instance.base = m_base;
    instance.index = index;           // Index of this instance
    instance.active = false;          // Used to determine if active
    instance.activecounter = 0;       // Counter used to flip active bit
    instance.activebase = 2. * pow(m_base,(index-1));

    // Storage for values of y currently used to update u
    instance.stage_y    = ComplexTripleArray( m_nQuadPts );
    instance.cstash_y   = ComplexTripleArray( m_nQuadPts );
    instance.csandbox_y = ComplexTripleArray( m_nQuadPts );
    instance.fstash_y   = ComplexTripleArray( m_nQuadPts );
    instance.fsandbox_y = ComplexTripleArray( m_nQuadPts );

    for( unsigned int q=0; q<m_nQuadPts; ++q )
    {
        instance.stage_y[q]    = ComplexDoubleArray( m_nvars );
        instance.cstash_y[q]   = ComplexDoubleArray( m_nvars );
        instance.csandbox_y[q] = ComplexDoubleArray( m_nvars );
        instance.fstash_y[q]   = ComplexDoubleArray( m_nvars );
        instance.fsandbox_y[q] = ComplexDoubleArray( m_nvars );

        for( unsigned int i=0; i<m_nvars; ++i )
        {
            instance.stage_y[q][i]    = ComplexSingleArray( m_npoints, 0.0 );
            instance.cstash_y[q][i]   = ComplexSingleArray( m_npoints, 0.0 );
            instance.csandbox_y[q][i] = ComplexSingleArray( m_npoints, 0.0 );
            instance.fstash_y[q][i]   = ComplexSingleArray( m_npoints, 0.0 );
            instance.fsandbox_y[q][i] = ComplexSingleArray( m_npoints, 0.0 );
        }
    }

    // Major storage for auxilliary ODE solutions.
    instance.stage_ind = std::pair<int, int>(0, 0);  // Time-step
                                                     // counters
                                                     // indicating to
                                                     // which interval
                                                     // ymain is
                                                     // associated

    // Staging allocation
    instance.stage_active = false;
    instance.stage_ccounter = 0;
    instance.stage_cbase = pow(m_base, index-1); // This base is halved
                                                 // after the first cycle
    instance.stage_fcounter = 0;
    instance.stage_fbase = pow(m_base, index);   // This base is halved
                                                 // after the first cycle

    // Ceiling stash allocation
    instance.cstash_counter = 0;               // Counter used to
                                               // determine when to
                                               // stash

    instance.cstash_base = pow(m_base, index-1);     // base for counter
    instance.cstash_ind = std::pair<int, int>(0, 0); // ind(1) is never
                                                     // used: it always
                                                     // matches
                                                     // main.ind(1)

    // Ceiling sandbox allocation
    instance.csandbox_active = false; // Flag to determine when stash 2
                                      // is utilized
    instance.csandbox_counter = 0;
    instance.csandbox_ind = std::pair<int, int>(0, 0);

    // Floor stash
    instance.fstash_base = 2*pow(m_base, index);
    instance.fstash_ind = std::pair<int, int>(0, 0);

    // Floor sandbox
    instance.fsandbox_active = false;
    instance.fsandbox_activebase = pow(m_base, index);
    instance.fsandbox_stashincrement = (m_base-1) * pow(m_base, index-1);
    instance.fsandbox_ind = std::pair<int, int>(0, 0);

    // Defining parameters of the Talbot contour quadrature rule
    NekDouble sigma = 0;
    NekDouble mu0 = 8;
    NekDouble nu = 0.6;
    NekDouble Tl =
      m_deltaT * (2.*pow(m_base, index) - 1. - pow(m_base, index-1));
    NekDouble mu = mu0 / Tl;

    // Talbot quadrature rule
    talbot_quadrature(m_nQuadPts, mu, nu, sigma, instance.z, instance.w);

    // With sigma == 0, the dependence of z and w on index is just a
    // multiplicative scaling factor (mu). So technically we'll only
    // need one instance of this N-point rule and can scale it
    // accordingly inside each integral_class instance. Not sure if
    // this optimization is worth it. Cumulative memory savings would
    // only be about 4*N*Lmax floats.

    // Below: precomputation for time integration of auxiliary
    // variables.  Everything below here is independent of the
    // instance index index. Therefore, we could actually just
    // generate and store one copy of this stuff and use it
    // everywhere.

    // std::cout << index << "  " << std::endl;

    instance.E = Array<OneD, std::complex<NekDouble> >(m_nQuadPts, 0.0);

    for( unsigned int q=0; q<m_nQuadPts; ++q )
    {
        instance.E[q] = exp(instance.z[q] * m_deltaT);

        // std::cout << instance.E[i] << "  ";
    }
    // std::cout << std::endl;

    instance.Eh =
        Array<TwoD, std::complex<NekDouble> >(m_order+1, m_nQuadPts, 0.0);

    for( unsigned int q=0; q<m_nQuadPts; ++q )
    {
      instance.Eh[0][q] =
        1. / instance.z[q] * (exp(instance.z[q] * m_deltaT) - 1.0);
    }

    for( unsigned int m=1; m<m_order+1; ++m )
    {
        for( unsigned int q=0; q<m_nQuadPts; ++q )
        {
            instance.Eh[m][q] = -1./instance.z[q] +
              NekDouble(m)/(instance.z[q] * m_deltaT) * instance.Eh[m-1][q];
        }
    }

    // std::cout << "instance.Eh " << std::endl;

    // for( unsigned int m=0; m<m_order+1; ++m )
    // {
    //   for( unsigned int q=0; q<m_nQuadPts/4; ++q )
    //   {
    //     std::cout << instance.Eh[m][q] << "  ";
    //   }
    //   std::cout << std::endl;
    // }
    // std::cout << std::endl;


    // As array one for each order.
    Array<OneD, Array<TwoD, NekDouble> > &As = instance.As;

    As = Array<OneD, Array<TwoD, NekDouble> >(m_order+2);

    for( unsigned int m=1; m<=m_order+1; ++m )
    {
        As[m] = Array<TwoD, NekDouble>(m,  m, 0.0);

        switch( m )
        {
          case 1:
            As[m][0][0] = 1.;
            break;

          case 2:
            As[m][0][0] =  0.;      As[m][0][1] =  1.;
            As[m][1][0] =  1.;      As[m][1][1] = -1.;
            break;

          case 3:
            As[m][0][0] =  0.;     As[m][0][1] =  1.;  As[m][0][2] =  0;
            As[m][1][0] =  1./2.;  As[m][1][1] =  0.;  As[m][1][2] = -1./2.;
            As[m][2][0] =  1./2.;  As[m][2][1] = -1.;  As[m][2][2] =  1./2.;
            break;

          case 4:
            As[m][0][0] =  0.;      As[m][0][1] =  1.;
            As[m][0][2] =  0.;      As[m][0][3] =  0.;

            As[m][1][0] =  1./3.;   As[m][1][1] = 1./2.;
            As[m][1][2] = -1.;      As[m][1][3] = 1./6.;

            As[m][2][0] =  1./2.;   As[m][2][1] = -1.;
            As[m][2][2] =  1./2.;   As[m][2][3] =  0.;

            As[m][3][0] =  1./6.;   As[m][3][1] = -1./2.;
            As[m][3][2] =  1./2.;   As[m][3][3] = -1./6.;
            break;

          case 5:
            As[m][0][0] =  0.;      As[m][0][1] =  1.;
            As[m][0][2] =  0.;      As[m][0][3] =  0.;    As[m][0][4] =  0.;

            As[m][1][0] =  1./4.;   As[m][1][1] =  5./6.;
            As[m][1][2] = -3./2.;   As[m][1][3] =  1./2.; As[m][1][4] = -1./12.;

            As[m][2][0] =  11./24.; As[m][2][1] = -5./6.;
            As[m][2][2] =  1./4.;   As[m][2][3] =  1./6.; As[m][2][4] = -1./24.;

            As[m][3][0] =  1./4.;   As[m][3][1] = -5./6.;
            As[m][3][2] =  1.;      As[m][3][3] = -1./2.; As[m][3][4] =  1./12.;

            As[m][4][0] =  1./24.;  As[m][4][1] = -1./6.;
            As[m][4][2] =  1./4.;   As[m][4][3] = -1./6.; As[m][4][4] =  1./24.;
            break;

            // The default is a general formula, but the matrix inversion
            // involved is ill-conditioned, so the special cases below are
            // epxlicitly given to combat roundoff error in the most-used
            // scenarios.
        default:
          ASSERTL1(false, "No matrix inverse.");

          // Ainv = zeros(counter);
          // Ainv(1,:) = 1;
          // Ainv(2,1) = 1;

          // for( unsigned int j = 3; j<=counter; ++j )
          // {
          //          Ainv(j,:) = pow((j-2)., (0:counter));
          //          Ainv(j,:) = Ainv(j,:). * pow((-1)., (0:counter));
          // }

          // As[m] = inv(Ainv);
          break;
        }
        // std::cout << m << std::endl;
        // for( unsigned int i=0; i<m; ++i ) {
        //   for( unsigned int j=0; j<m; ++j ) {
        //  std::cout << As[m][i][j] << "  ";
        //   }
        //   std::cout << std::endl;
        // }

        // std::cout << std::endl;
    }

    instance.AtEh  =
      Array<TwoD, std::complex<NekDouble> >(m_order+1, m_nQuadPts, 0.0);
    // instance.AptEh =
    //   Array<TwoD, std::complex<NekDouble> >(m_order,   m_nQuadPts, 0.0);

    // std::cout << "Instance.As" << std::endl;
    // for( unsigned int k=0; k<m_order+1; ++k )
    // {
    //   for( unsigned int j=0; j<m_order+1; ++j )
    //   {
    //     std::cout << As[m_order+1][k][j] << "  " ;
    //   }

    //   std::cout << std::endl;
    // }
    // std::cout << std::endl;

    // std::cout << "Instance.Eh" << std::endl;
    // for( unsigned int k=0; k<m_order+1; ++k )
    // {
    //   for( unsigned int q=0; q<m_nQuadPts/4; ++q )
    //   {
    //     std::cout << instance.Eh[k][q] << "  " ;
    //   }

    //   std::cout << std::endl;
    // }
    // std::cout << std::endl;

    for( unsigned int m=0; m<m_order+1; ++m )
    {
        for( unsigned int q=0; q<m_nQuadPts; ++q )
        {
            for( unsigned int i=0; i<m_order+1; ++i )
            {
              instance.AtEh[m][q] += As[m_order+1][m][i] * instance.Eh[i][q];
            }
        }

      // for( unsigned int j=0; j<m_order; ++j )
      // {
      //        for( unsigned int i=0; i<m_order; ++i )
      //        {
      //          instance.AptEh[j][q] += As[m_order][j][i] * instance.Eh[i][q];
      //        }
      // }
    }

    // std::cout << "instance.AtEh" << std::endl;
    // for( unsigned int j=0; j<m_order+1; ++j )
    // {
    //   for( unsigned int q=0; q<m_nQuadPts/4; ++q )
    //   {
    //  std::cout << instance.AtEh[j][q] << "  ";
    //   }

    //   std::cout << std::endl;
    // }

    // std::cout << std::endl;

    // for( unsigned int j=0; j<m_order+1; ++j )
    // {
    //   for( unsigned int q=0; q<m_nQuadPts/4; ++q )
    //   {
    //  std::cout << instance.AptEh[j][q] << "  ";
    //   }
    // }

    // std::cout << std::endl;
}


// Performs rearrangement of staging/stashing for current time
//
// (1) activates ceiling staging
// (2) moves ceiling stash ---> stage
// (3) moves floor stash --> stage (+ updates all ceiling data)
void FractionalInTimeIntegrationScheme::
update_stage(const unsigned int timeStep,
             Instance &instance)
{
    // Counter to flip active bit
    if( !instance.active )
    {
        instance.active = (timeStep % instance.activebase == 0);
    }
    // std::cout << timeStep << "  "
    //        << instance.base << "  "
    //        << instance.index << "  "
    //        << instance.active << "  "
    //        << instance.activecounter << "  "
    //        << instance.activebase << "  "
    //        << std::endl;

    // Determine if staging is necessary
    if( instance.active )
    {
      // std::cout << timeStep << "     "
      //      << instance.base << "  "
      //      << instance.index << "  "
      //      << instance.active << "  "
      //      << instance.activecounter << "  "
      //      << instance.activebase << "  "
      //      << std::endl;

        // Floor staging superscedes ceiling staging
        instance.stage_ind = instance.fstash_ind;
        // instance.stage_y   = instance.fstash_y;
        if( timeStep % instance.fstash_base == 0 )
        {
            for( unsigned int q=0; q<m_nQuadPts; ++q )
            {
                for( unsigned int i=0; i<m_nvars; ++i )
                {
                    for( unsigned int j=0; j<m_npoints; ++j )
                    {
                        instance.stage_y[q][i][j] = instance.fstash_y[q][i][j];
                    }
                }
            }

            // std::cout << "Floor staging superscedes ceiling staging " << instance.index << std::endl;
            // for( unsigned int q=0; q<m_nQuadPts/4; ++q )
            //        {
            //          for( unsigned int i=0; i<m_nvars; ++i )
            //            {
            //              for( unsigned int j=0; j<m_npoints; ++j )
            //                {
            //                  std::cout << instance.stage_y[q][i][j] << "  ";
            //                }
            //              std::cout << std::endl;
            //            }
            //          std::cout << std::endl;
            //        }
            // std::cout << std::endl;

            instance.csandbox_ind = instance.fsandbox_ind;
            // instance.csandbox_y   = instance.fsandbox_y;
            for( unsigned int q=0; q<m_nQuadPts; ++q )
            {
                for( unsigned int i=0; i<m_nvars; ++i )
                {
                    for( unsigned int j=0; j<m_npoints; ++j )
                    {
                        instance.csandbox_y[q][i][j] = instance.fsandbox_y[q][i][j];
                    }
                }
            }

            // After floor staging happens once, new base is base^index
            instance.fstash_base = pow(instance.base, instance.index);

            // Restart floor sandbox
            // instance.fsandbox_y = zeros(size(instance.fsandbox_y));
            for( unsigned int q=0; q<m_nQuadPts; ++q )
            {
                for( unsigned int i=0; i<m_nvars; ++i )
                {
                    for( unsigned int j=0; j<m_npoints; ++j )
                    {
                        instance.fsandbox_y[q][i][j] = 0;
                    }
                }
            }

            instance.fsandbox_ind = std::pair<int, int>(0,0);
            instance.fsandbox_active = false;
        }
        // Check for ceiling staging
        else if( timeStep % instance.stage_cbase == 0 )
        {
            instance.stage_ind = instance.cstash_ind;
            // instance.stage_y   = instance.cstash_y;
            for( unsigned int q=0; q<m_nQuadPts; ++q )
            {
                for( unsigned int i=0; i<m_nvars; ++i )
                {
                    for( unsigned int j=0; j<m_npoints; ++j )
                    {
                        instance.stage_y[q][i][j] = instance.cstash_y[q][i][j];
                    }
                }
            }

            // std::cout << "Check for ceiling staging " << instance.index << std::endl;
            // for( unsigned int q=0; q<m_nQuadPts/4; ++q )
            //        {
            //          for( unsigned int i=0; i<m_nvars; ++i )
            //            {
            //              for( unsigned int j=0; j<m_npoints; ++j )
            //                {
            //                  std::cout << instance.stage_y[q][i][j] << "  ";
            //                }
            //              std::cout << std::endl;
            //            }
            //          std::cout << std::endl;
            //        }
            // std::cout << std::endl;
        }
    }
}

// Approximates the integral
//
//   \int_{(m-1) h}^{m h} k(m*h -s) f(u, s) dx{s}
//
// using a time-stepping scheme of a particular order. Here, k depends on alph,
// the derivative order.
void FractionalInTimeIntegrationScheme::
final_increment(const unsigned int timeStep,
                const TimeIntegrationSchemeOperators &op)
{
    for( unsigned int m=0; m<m_order; ++m )
    {
        op.DoOdeRhs(m_u[m], m_F[m], m_deltaT * (timeStep-m));
    }

    for( unsigned int i=0; i<m_nvars; ++i )
    {
        for( unsigned int j=0; j<m_npoints; ++j )
        {
            m_uNext[i][j] = 0;

            for( unsigned int q=0; q<m_order; ++q )
            {
                m_uNext[i][j] += m_F[q][i][j] * m_AhattJ[q];
            }
            // std::cout << m_uNext[i][j];
        }
        // std::cout << std::endl;
    }
    // std::cout << std::endl;
}


void FractionalInTimeIntegrationScheme::
integral_contribution(const unsigned int timeStep,
                      const unsigned int tauml,
                      const Instance &instance)
{
    // Assume y has already been updated to time level m
    for( unsigned int q=0; q<m_nQuadPts; ++q )
    {
        m_expFactor[q] =
            exp(instance.z[q] * m_deltaT * NekDouble(timeStep - tauml)) *
            pow(instance.z[q], -m_alpha) * instance.w[q];
    }

    // std::cout << "integral_contribution instance.stage_y " << instance.index << std::endl;
    // for( unsigned int i=0; i<m_nvars; ++i )
    //   {
    //     for( unsigned int q=0; q<m_nQuadPts/4; ++q )
    //    {
    //      std::cout << instance.stage_y[i][q] << "  ";
    //    }
    //     std::cout << std::endl;
    //   }
    // std::cout << std::endl;

    for( unsigned int i=0; i<m_nvars; ++i )
    {
        for( unsigned int j=0; j<m_npoints; ++j )
        {
            m_uInt[i][j] = 0;

            for( unsigned int q=0; q<m_nQuadPts; ++q )
            {
                m_uInt[i][j] += instance.stage_y[q][i][j] * m_expFactor[q];
            }

            if( m_uInt[i][j].real() < 1e8 )
            {
                m_uInt[i][j] = m_uInt[i][j].real();
            }
        }
    }

    // std::cout << "Instance " << instance.index << std::endl;
    // for( unsigned int i=0; i<m_nvars; ++i )
    //   {
    //     for( unsigned int j=0; j<m_npoints; ++j )
    //    {
    //        std::cout << I[i] << "  ";
    //    }
    //     std::cout << std::endl;
    //   }
    // std::cout << std::endl;
}


void FractionalInTimeIntegrationScheme::
time_advance(const unsigned int timeStep,
             const Instance &instance,
             const TimeIntegrationSchemeOperators &op,
             ComplexTripleArray &y)
{
    // Solution to y' = z*y + f(u), using an exponential integrator with
    // implicit order (m_order + 1) interpolation of the f(u) term.

    // std::cout << instance.index << std::endl;

    // Try automated high-order method:
    if( timeStep < m_order )
    {
        // Not enough history. For now, demote to lower-order method.
        // TODO: use multi-stage method
        for( unsigned int m=0; m<timeStep+1; ++m )
        {
            op.DoOdeRhs(m_u[m], m_F[m], m_deltaT * (timeStep-m));
        }

        // std::cout << "F" << std::endl;
        // for( unsigned int m=0; m<timeStep+1; ++m )
        //   {
        //  for( unsigned int i=0; i<m_nvars; ++i )
        //    {
        //      for( unsigned int j=0; j<m_npoints; ++j )
        //        {
        //          std::cout << m_F[m][i][j] << "  ";
        //        }
        //      std::cout << std::endl;
        //    }
        //  std::cout << std::endl;
        //   }
        // std::cout << std::endl;

        // y = y * instance.E + F * instance.As{timeStep+1}.' * instance.Eh(1:(timeStep+1),:);

        // instance.As{timeStep+1}.' * instance.Eh(1:(timeStep+1),:);
        Array<TwoD, std::complex<NekDouble> > tmp =
          Array<TwoD, std::complex<NekDouble> >(timeStep+1, m_nQuadPts, 0.0);

        for( unsigned int m=0; m<timeStep+1; ++m )
        {
            for( unsigned int q=0; q<m_nQuadPts; ++q )
            {
                for( unsigned int i=0; i<timeStep+1; ++i )
                {
                    tmp[m][q] +=
                        instance.As[timeStep+1][m][i] * instance.Eh[i][q];
                }
            }
        }

        // std::cout << "tmp" << std::endl;
        // for( unsigned int k=0; k<timeStep+1; ++k )
        //   {
        //  for( unsigned int q=0; q<m_nQuadPts/4; ++q )
        //    {
        //      std::cout << tmp[k][q] << "  ";
        //    }
        //  std::cout << std::endl;
        //   }
        // std::cout << std::endl;

        // y = y * instance.E + F * instance.As{timeStep+1}.' * instance.Eh(1:(timeStep+1),:);
        for( unsigned int q=0; q<m_nQuadPts; ++q )
        {
            for( unsigned int i=0; i<m_nvars; ++i )
            {
                for( unsigned int j=0; j<m_npoints; ++j )
                {
                    // y * instance.E
                    y[q][i][j] *= instance.E[q];

                    for( unsigned int m=0; m<timeStep+1; ++m )
                    {
                        // F * instance.As{timeStep+1}.' * instance.Eh(1:(timeStep+1),:);
                        y[q][i][j] += m_F[m][i][j] * tmp[m][q];
                    }
                }
            }
        }
    }
    else
    {
        for( unsigned int m=0; m<m_order+1; ++m )
        {
            op.DoOdeRhs(m_u[m], m_F[m], m_deltaT * (timeStep-m));
        }

        // std::cout << "F" << std::endl;
        // for( unsigned int m=0; m<m_order+1; ++m )
        //   {
        //  for( unsigned int i=0; i<m_nvars; ++i )
        //    {
        //      for( unsigned int j=0; j<m_npoints; ++j )
        //        {
        //          std::cout << m_F[m][i][j] << "  ";
        //        }
        //      std::cout << std::endl;
        //    }
        //  std::cout << std::endl;
        //   }
        // std::cout << std::endl;

        // y = y * instance.E + F * instance.AtEh;
        for( unsigned int q=0; q<m_nQuadPts; ++q )
        {
            for( unsigned int i=0; i<m_nvars; ++i )
            {
                for( unsigned int j=0; j<m_npoints; ++j )
                {
                    // y * instance.E
                    y[q][i][j] *= instance.E[q];

                    for( unsigned int m=0; m<m_order+1; ++m )
                    {
                        // F * instance.AtEh
                        y[q][i][j] += m_F[m][i][j] * instance.AtEh[m][q];
                    }
                }
            }
        }
    }

    // std::cout << "y" << std::endl;
    // for( unsigned int q=0; q<m_nQuadPts/4; ++q )
    //   {
    //    for( unsigned int i=0; i<m_nvars; ++i )
    //      {
    //        for( unsigned int j=0; j<m_npoints; ++j )
    //          {
    //            std::cout << y[q][i][j] << "  ";
    //          }
    //        std::cout << std::endl;
    //      }
    //    std::cout << std::endl;
    //   }
    // std::cout << std::endl;
}

// Updates sandboxes to current time
// (1) advances ceiling sandbox
// (2) moves ceiling sandbox ---> stash
// (3) activates floor sandboxing
// (4) advances floor sandbox
// (5) moves floor sandbox ---> stash

void FractionalInTimeIntegrationScheme::
advance_sandbox(const unsigned int timeStep,
                const TimeIntegrationSchemeOperators &op,
                Instance &instance)
{
    // (1)
    // update(instance.csandbox.y)
    time_advance(timeStep, instance, op, instance.csandbox_y);
    instance.csandbox_ind.second = timeStep;

    // std::cout << "instance.csandbox_y " << instance.index << std::endl;
    // for( unsigned int q=0; q<m_nQuadPts/4; ++q )
    //   {
    //     for( unsigned int i=0; i<m_nvars; ++i )
    //    {
    //      for( unsigned int j=0; j<m_npoints; ++j )
    //        {
    //          std::cout << instance.csandbox_y[q][i][j] << "  ";
    //        }
    //      std::cout << std::endl;
    //    }
    //     std::cout << std::endl;
    //   }
    // std::cout << std::endl;

    // (2)
    // Determine if ceiling stashing is necessary
    instance.cstash_counter = modincrement(instance.cstash_counter,
                                           instance.cstash_base);

    // std::cout << instance.index << "  " << instance.cstash_counter << std::endl;

    if( timeStep % instance.cstash_base == 0 )
    {
        // std::cout << "Then need to stash" << std::endl;
        // Then need to stash
        instance.cstash_ind = instance.csandbox_ind;
        // instance.cstash_y   = instance.csandbox_y;
        for( unsigned int q=0; q<m_nQuadPts; ++q )
        {
            for( unsigned int i=0; i<m_nvars; ++i )
            {
                for( unsigned int j=0; j<m_npoints; ++j )
                {
                    instance.cstash_y[q][i][j] = instance.csandbox_y[q][i][j];
                }
            }
        }
    }

    if( instance.fsandbox_active )
    {
        // std::cout << "instance.fsandbox_active" << std::endl;
        // (4)
        time_advance(timeStep, instance, op, instance.fsandbox_y);

        // std::cout << "instance.fsandbox_y " << instance.index << std::endl;
        // for( unsigned int q=0; q<m_nQuadPts/4; ++q )
        //   {
        //  for( unsigned int i=0; i<m_nvars; ++i )
        //    {
        //      for( unsigned int j=0; j<m_npoints; ++j )
        //        {
        //          std::cout << instance.fsandbox_y[q][i][j] << "  ";
        //        }
        //      std::cout << std::endl;
        //    }
        //  std::cout << std::endl;
        //   }
        // std::cout << std::endl;

        instance.fsandbox_ind.second = timeStep;

        // (5) Move floor sandbox to stash
        if( (instance.fsandbox_ind.second - instance.fsandbox_ind.first) %
            instance.fsandbox_stashincrement == 0 )
        {
            // std::cout << "Move floor sandbox to stash" << std::endl;
            instance.fstash_ind = instance.fsandbox_ind;
            // instance.fstash_y   = instance.fsandbox_y;
            for( unsigned int q=0; q<m_nQuadPts; ++q )
            {
                for( unsigned int i=0; i<m_nvars; ++i )
                {
                    for( unsigned int j=0; j<m_npoints; ++j )
                    {
                        instance.fstash_y[q][i][j] =
                            instance.fsandbox_y[q][i][j];
                    }
                }
            }
        }
    }
    else // Determine if advancing floor sandbox is necessary at next time
    {
        // (3)
        if( timeStep % instance.fsandbox_activebase == 0 )
        {
            // std::cout << "Determine if advancing floor sandbox is necessary at next time" << std::endl;
            instance.fsandbox_active = true;
            instance.fsandbox_ind = std::pair<int,int>(timeStep, timeStep);
        }
    }
}

} // end namespace LibUtilities
} // end namespace NekTar
