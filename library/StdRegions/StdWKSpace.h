#ifndef STDWKSPACE_H
#define STDWKSPACE_H

#include <LibUtilities/ErrorUtil.hpp>

#include <StdRegions/StdRegions.hpp>

namespace Nektar
{
namespace StdRegions
{

  class WorkSpace
  {
  public:
    WorkSpace(){
      _availspace=0;
    }

    ~WorkSpace(){
    }

    /// Get a temporary workspace of size SIZE_WSTHREADS
    inline double * get_workspace(){
      if(_availspace < NUM_WSTHREADS*NUM_STDWKSPC){
    _availspace++;
    return(&_workspace[_availspace-1][0]);
      }
      else
    ErrorUtil::Error(fatal,"WorkSpace", "Workspace availability exceeded");
      return NULL;
    }

    /// Free a temporary workspace
    inline void free_workspace(){
      _availspace--;
      if(_availspace < 0)
    ErrorUtil::Error(fatal,"WorkSpace", "Invalid workspace free");
    }

  private:
    double _workspace[NUM_WSTHREADS*NUM_STDWKSPC][SIZE_WSTHREADS];
    int _availspace;

  };

  class WorkSpaceHandle{
  private:
    double * _ptr;
    static WorkSpace _WS; ///< Workspace class used by StdExpansions

  public:
    WorkSpaceHandle(){
      _ptr = _WS.get_workspace();
    }

    ~WorkSpaceHandle(){
      _ptr = NULL;
      _WS.free_workspace();
    }

    inline double * get_workspace(){
      return(_ptr);
    }

    inline double operator[](const int i) const{
      ASSERTL2(i>=0,"A: [] index less than zero");
      ASSERTL2(i<SIZE_WSTHREADS,
          "A: [] index greater than SIZE\_WSTHREADS");
      return _ptr[i];
    }

    inline double& operator[](const int i){
      ASSERTL2(i>=0,"B: [] index less than zero");
      ASSERTL2(i<SIZE_WSTHREADS,
           "B: [] index greater than SIZE\_WSTHREADS");
      return _ptr[i];
    }
  };

}//end namespace
}//end namespace

#endif  //STDWKSPACE_H
