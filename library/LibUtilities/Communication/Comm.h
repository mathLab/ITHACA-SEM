///////////////////////////////////////////////////////////////////////////////
//
// File Comm.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: Base communication class
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_UTILITIES_COMM_H
#define NEKTAR_LIB_UTILITIES_COMM_H

#include <vector>

#include <boost/enable_shared_from_this.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
//namespace Nektar { template <typename Dim, typename DataType> class Array; }
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Communication/CommDataType.h>

namespace Nektar
{
    namespace LibUtilities
    {
        // Forward declarations
        class Comm;

        /// Pointer to a Communicator object.
        typedef boost::shared_ptr<Comm> CommSharedPtr;

        /// Datatype of the NekFactory used to instantiate classes derived from
        /// the EquationSystem class.
        typedef LibUtilities::NekFactory< std::string, Comm, int, char** > CommFactory;

        LIB_UTILITIES_EXPORT CommFactory& GetCommFactory();

		/// Type of operation to perform in AllReduce.
        enum ReduceOperator
        {
            ReduceSum,
            ReduceMax,
            ReduceMin
        };

        /// Base communications class
        class Comm: public boost::enable_shared_from_this<Comm>
        {
            public:
                LIB_UTILITIES_EXPORT Comm(int narg, char* arg[]);
                LIB_UTILITIES_EXPORT virtual ~Comm();

                LIB_UTILITIES_EXPORT inline void Finalise();

                /// Returns number of processes
                LIB_UTILITIES_EXPORT inline int GetSize();
                LIB_UTILITIES_EXPORT inline int GetRank();
                LIB_UTILITIES_EXPORT inline const std::string& GetType() const;

                /// Block execution until all processes reach this point
                LIB_UTILITIES_EXPORT inline void Block();

                template<class T>
                void Send(int pProc, const T& pData);
                template<class T>
                void Recv(int pProc, T& pData);
                template<class T>
                void SendRecv(int pSendProc, const T& pSendData,
                        int pRecvProc, T& pRecvData);
                template<class T>
                void SendRecvReplace(int pSendProc, int pRecvProc, T& pData);

                template<class T>
                void AllReduce(T& pData, enum ReduceOperator pOp);

                LIB_UTILITIES_EXPORT inline void AlltoAll(Array<OneD, NekDouble>& pSendData,
                                                          Array<OneD, NekDouble>& pRecvData);
                LIB_UTILITIES_EXPORT inline void AlltoAll(Array<OneD, int>& pSendData,
                                                          Array<OneD, int>& pRecvData);
                LIB_UTILITIES_EXPORT inline void AlltoAllv(Array<OneD, NekDouble>& pSendData,
                                                           Array<OneD, int>& pSendDataSizeMap,
                                                           Array<OneD, int>& pSendDataOffsetMap,
                                                           Array<OneD, NekDouble>& pRecvData,
                                                           Array<OneD, int>& pRecvDataSizeMap,
                                                           Array<OneD, int>& pRecvDataOffsetMap);
                LIB_UTILITIES_EXPORT inline void AlltoAllv(Array<OneD, int>& pSendData,
                                                           Array<OneD, int>& pSendDataSizeMap,
                                                           Array<OneD, int>& pSendDataOffsetMap,
                                                           Array<OneD, int>& pRecvData,
                                                           Array<OneD, int>& pRecvDataSizeMap,
                                                           Array<OneD, int>& pRecvDataOffsetMap);
                
                template<class T>
                void Bcast(T& data, int rootProc);

                LIB_UTILITIES_EXPORT inline void Exscan(const Array<OneD, unsigned long long>& pData, const enum ReduceOperator pOp, Array<OneD, unsigned long long>& ans);

                LIB_UTILITIES_EXPORT inline Array<OneD, unsigned long long> Gather(const int rootProc, const Array<OneD, unsigned long long>& val);
                LIB_UTILITIES_EXPORT inline Array<OneD, unsigned long long> Scatter(const int rootProc, const Array<OneD, unsigned long long>& pData);

                LIB_UTILITIES_EXPORT inline CommSharedPtr CommCreateIf(int flag);

                LIB_UTILITIES_EXPORT inline void SplitComm(int pRows, int pColumns);
                LIB_UTILITIES_EXPORT inline CommSharedPtr GetRowComm();
                LIB_UTILITIES_EXPORT inline CommSharedPtr GetColumnComm();


                LIB_UTILITIES_EXPORT inline bool TreatAsRankZero(void);
                LIB_UTILITIES_EXPORT inline bool RemoveExistingFiles(void);

            protected:
                int m_size;                 ///< Number of processes
                std::string m_type;         ///< Type of communication
                CommSharedPtr m_commRow;    ///< Row communicator
                CommSharedPtr m_commColumn; ///< Column communicator

                Comm();

                virtual void v_Finalise() = 0;
                virtual int  v_GetRank() = 0;
                virtual void v_Block() = 0;
                virtual void v_Send(const void* buf, int count, CommDataType dt, int dest) = 0;
                virtual void v_Recv(void* buf, int count, CommDataType dt, int source) = 0;
                virtual void v_Sendrecv(const void *sendbuf, int sendcount, CommDataType sendtype, int dest,
                        void *recvbuf, int recvcount, CommDataType recvtype, int source) = 0;
				virtual void v_SendRecvReplace(void* buf, int count, CommDataType dt,
				        int pSendProc, int pRecvProc) = 0;
				virtual void v_AllReduce(void* buf, int count, CommDataType dt, enum ReduceOperator pOp) = 0;

				virtual void v_AlltoAll(Array<OneD, NekDouble>& pSendData,
										Array<OneD, NekDouble>& pRecvData) = 0;
                virtual void v_AlltoAll(Array<OneD, int>& pSendData,
										Array<OneD, int>& pRecvData) = 0;
			    virtual void v_AlltoAllv(Array<OneD, NekDouble>& pSendData,
										Array<OneD, int>& pSendDataSizeMap,
										Array<OneD, int>& pSendDataOffsetMap,
										Array<OneD, NekDouble>& pRecvData,
										Array<OneD, int>& pRecvDataSizeMap,
										Array<OneD, int>& pRecvDataOffsetMap) = 0;
				virtual void v_AlltoAllv(Array<OneD, int>& pSendData,
										Array<OneD, int>& pSendDataSizeMap,
										Array<OneD, int>& pSendDataOffsetMap,
										Array<OneD, int>& pRecvData,
										Array<OneD, int>& pRecvDataSizeMap,
										Array<OneD, int>& pRecvDataOffsetMap) = 0;
				virtual void v_Bcast(void* buffer, int count, CommDataType dt, int root) = 0;
                virtual void v_Exscan(const Array<OneD, unsigned long long>& pData, const enum ReduceOperator pOp, Array<OneD, unsigned long long>& ans) = 0;

                virtual Array<OneD, unsigned long long> v_Gather(const int rootProc, const Array<OneD, unsigned long long>& val) = 0;
                virtual Array<OneD, unsigned long long> v_Scatter(const int rootProc, const Array<OneD, unsigned long long>& pData) = 0;

                virtual CommSharedPtr v_CommCreateIf(int flag) = 0;
                virtual void v_SplitComm(int pRows, int pColumns) = 0;
                virtual bool v_TreatAsRankZero(void) = 0;
                LIB_UTILITIES_EXPORT virtual bool v_RemoveExistingFiles(void);
        };


        /**
         *
         */
        inline void Comm::Finalise()
        {
            v_Finalise();
        }

        /**
         *
         */
        inline int Comm::GetSize()
        {
            return m_size;
        }

        /**
         *
         */
        inline int Comm::GetRank()
        {
            return v_GetRank();
        }

        /**
         *
         */
        inline const std::string& Comm::GetType() const
        {
            return m_type;
        }

        /**
         *
         */
        inline void Comm::Block()
        {
            v_Block();
        }

        template<class T>
        void Comm::Send(int pProc, const T& pData)
        {
                v_Send(CommDataTypeTraits<T>::GetPointer(pData),
                        CommDataTypeTraits<T>::GetCount(pData),
                        CommDataTypeTraits<T>::GetDataType(),
                        pProc);
        }

        template<class T>
        void Comm::Recv(int pProc, T& pData)
        {
                v_Recv(CommDataTypeTraits<T>::GetPointer(pData),
                        CommDataTypeTraits<T>::GetCount(pData),
                        CommDataTypeTraits<T>::GetDataType(),
                        pProc);
        }

        /**
         *
         */
        template<class T>
        void Comm::SendRecv(int pSendProc,
                             const T& pSendData,
                             int pRecvProc,
                             T& pRecvData)
        {
            v_SendRecv(CommDataTypeTraits<T>::GetPointer(pSendData),
                    CommDataTypeTraits<T>::GetCount(pSendData).
                    CommDataTypeTraits<T>::GetDataType(),
                    pSendProc,
                    CommDataTypeTraits<T>::GetPointer(pRecvData),
                    CommDataTypeTraits<T>::GetCount(pRecvData),
                    CommDataTypeTraits<T>::GetDataType(),
                    pRecvProc);
        }
		
		/**
         *
         */
        template <class T>
        void Comm::SendRecvReplace(int pSendProc,
                int pRecvProc,
                T& pData)
        {
            v_SendRecvReplace(CommDataTypeTraits<T>::GetPointer(pData),
                    CommDataTypeTraits<T>::GetCount(pData),
                    CommDataTypeTraits<T>::GetDataType(),
                    pSendProc, pRecvProc);
        }
		
        /**
         *
         */
        template <class T>
        void Comm::AllReduce(T& pData, enum ReduceOperator pOp)
        {
            v_AllReduce(CommDataTypeTraits<T>::GetPointer(pData),
                    CommDataTypeTraits<T>::GetCount(pData),
                    CommDataTypeTraits<T>::GetDataType(),
                    pOp);
        }

        /**
         *
         */
		inline void Comm::AlltoAll(Array<OneD, NekDouble>& pSendData,Array<OneD, NekDouble>& pRecvData)
		{
			v_AlltoAll(pSendData,pRecvData);
		}
		
		
		/**
         *
         */
		inline void Comm::AlltoAll(Array<OneD, int>& pSendData,Array<OneD, int>& pRecvData)
		{
			v_AlltoAll(pSendData,pRecvData);
		}
		
		
		/**
         *
         */
		inline void Comm::AlltoAllv(Array<OneD, NekDouble>& pSendData,
								 Array<OneD, int>& pSendDataSizeMap,
								 Array<OneD, int>& pSendDataOffsetMap,
								 Array<OneD, NekDouble>& pRecvData,
								 Array<OneD, int>& pRecvDataSizeMap,
								 Array<OneD, int>& pRecvDataOffsetMap)
		{
			v_AlltoAllv(pSendData,pSendDataSizeMap,pSendDataOffsetMap,pRecvData,pRecvDataSizeMap,pRecvDataOffsetMap);
		}
		
		/**
         *
         */
		inline void Comm::AlltoAllv(Array<OneD, int>& pSendData,
								 Array<OneD, int>& pSendDataSizeMap,
								 Array<OneD, int>& pSendDataOffsetMap,
								 Array<OneD, int>& pRecvData,
								 Array<OneD, int>& pRecvDataSizeMap,
								 Array<OneD, int>& pRecvDataOffsetMap)
		{
			v_AlltoAllv(pSendData,pSendDataSizeMap,pSendDataOffsetMap,pRecvData,pRecvDataSizeMap,pRecvDataOffsetMap);
		}

		template<class T>
		void Comm::Bcast(T& pData, int pRoot)
		{
		        v_Bcast(CommDataTypeTraits<T>::GetPointer(pData),
		                CommDataTypeTraits<T>::GetCount(pData),
		                CommDataTypeTraits<T>::GetDataType(),
		                pRoot);
		}

        inline void Comm::Exscan(const Array<OneD, unsigned long long>& pData, const enum ReduceOperator pOp, Array<OneD, unsigned long long>& ans)
        {
            v_Exscan(pData, pOp, ans);
        }
        /**
         * Concatenate all the input arrays, in rank order, onto the process with rank == rootProc
         */
        inline Array<OneD, unsigned long long> Comm::Gather(const int rootProc, const Array<OneD, unsigned long long>& val)
        {
            return v_Gather(rootProc, val);
        }
        /**
         * Scatter pData across ranks in chunks of len(pData)/num_ranks
         */
        inline Array<OneD, unsigned long long> Comm::Scatter(const int rootProc, const Array<OneD, unsigned long long>& pData)
        {
            return v_Scatter(rootProc, pData);
        }

		/**
		 * @brief If the flag is non-zero create a new communicator.
		 */
		inline CommSharedPtr Comm::CommCreateIf(int flag)
		{
		    return v_CommCreateIf(flag);
		}

        /**
         * @brief Splits this communicator into a grid of size pRows*pColumns
         * and creates row and column communicators. By default the communicator
         * is a single row.
         */
        inline void Comm::SplitComm(int pRows, int pColumns)
        {
            v_SplitComm(pRows, pColumns);
        }


        /**
         * @brief Retrieve the row communicator to which this process belongs.
         */
        inline CommSharedPtr Comm::GetRowComm()
        {
            if (!m_commRow.get())
            {
                return shared_from_this();
            }
            else
            {
                return m_commRow;
            }
        }


        /**
         * @brief Retrieve the column communicator to which this process
         * belongs.
         */
        inline CommSharedPtr Comm::GetColumnComm()
        {
            if (!m_commColumn.get())
            {
                return shared_from_this();
            }
            else
            {
                return m_commColumn;
            }
        }

        inline bool Comm::TreatAsRankZero(void)
        {
            return v_TreatAsRankZero();
        }

        inline bool Comm::RemoveExistingFiles(void)
        {
            return v_RemoveExistingFiles();
        }


    }
}

#endif
