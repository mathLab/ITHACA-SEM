#include <boost/algorithm/string.hpp>


#include <IncNavierStokesSolver/EquationSystems/CoupledLinearNS_sparse.h>
//include <IncNavierStokesSolver/EquationSystems/CoupledLinearNS.h>

#include <LibUtilities/BasicUtils/Timer.h>
#include <LocalRegions/MatrixKey.h>
#include <MultiRegions/GlobalLinSysDirectStaticCond.h>
#include <MultiRegions/ContField.h>

using namespace std;

namespace Nektar
{

    string CoupledLinearNS_sparse::className = SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction("CoupledLinearNS_sparse", CoupledLinearNS_sparse::create);


    /**
     *  @class CoupledLinearNS
     *
     * Set up expansion field for velocity and pressure, the local to
     * global mapping arrays and the basic memory definitions for
     * coupled matrix system
     */
    CoupledLinearNS_sparse::CoupledLinearNS_sparse(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
        : UnsteadySystem(pSession, pGraph),
          CoupledLinearNS(pSession, pGraph),
          m_zeroMode(false)
    {
    }



    void CoupledLinearNS_sparse::v_DoInitialise(void)
    {
		cout << "Hello from CoupledLinearNS_sparse::v_DoInitialise(void)" << endl;

    	load_session_parameters();
		load_snapshots();


	}

    void CoupledLinearNS_sparse::load_snapshots()
    {
    
    	cout << " Loading FOM snapshots from files ... " << endl;
    	
	// fill the fields snapshot_x_collection and snapshot_y_collection

	int nvelo = 2;
        Array<OneD, Array<OneD, NekDouble> > test_load_snapshot(nvelo); // for a 2D problem

        for(int i = 0; i < nvelo; ++i)
        {
            test_load_snapshot[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);  // number of phys points
        }
               
        std::vector<std::string> fieldStr;
        for(int i = 0; i < nvelo; ++i)
        {
           fieldStr.push_back(m_session->GetVariable(i));
        }

	snapshot_x_collection = Array<OneD, Array<OneD, NekDouble> > (number_of_snapshots);
	snapshot_y_collection = Array<OneD, Array<OneD, NekDouble> > (number_of_snapshots);

        for(int i = 0; i < number_of_snapshots; ++i)
        {
		// generate the correct string
		std::stringstream sstm;
		sstm << "TestSnap" << i+1;
		std::string result = sstm.str();
		const char* rr = result.c_str();

//	        EvaluateFunction(fieldStr, test_load_snapshot, result);
	        GetFunction(result)->Evaluate(fieldStr, test_load_snapshot);
	        
            //     now:   GetFunction( "Restart")->Evaluate(fieldStr,  Restart);
            //     PREV:    EvaluateFunction(fieldStr, Restart, "Restart");
	        
		snapshot_x_collection[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
		snapshot_y_collection[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
		for (int j=0; j < GetNpoints(); ++j)
		{
			snapshot_x_collection[i][j] = test_load_snapshot[0][j];
			snapshot_y_collection[i][j] = test_load_snapshot[1][j];
		}
        }    	

	
    	cout << " ... finished loading FOM snapshots from files" << endl;
    	
    	if (use_fine_grid_VV_and_load_ref && use_fine_grid_VV)
    	{
    		cout << "start loading Verification & Validation snapshots ... " << endl;
		snapshot_x_collection_VV = Array<OneD, Array<OneD, NekDouble> > (fine_grid_dir0);
		snapshot_y_collection_VV = Array<OneD, Array<OneD, NekDouble> > (fine_grid_dir0);    	
	        for(int i = 0; i < fine_grid_dir0; ++i)
	        {
			// generate the correct string
			std::stringstream sstm;
			sstm << "VV" << i+1;
			std::string result = sstm.str();
			const char* rr = result.c_str();
		        GetFunction(result)->Evaluate(fieldStr, test_load_snapshot);
		        
			snapshot_x_collection_VV[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
			snapshot_y_collection_VV[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
			for (int j=0; j < GetNpoints(); ++j)
			{
				snapshot_x_collection_VV[i][j] = test_load_snapshot[0][j];
				snapshot_y_collection_VV[i][j] = test_load_snapshot[1][j];
			}
	        }    	
    		cout << " ... finished loading Verification & Validation snapshots " << endl;
    	}
    }

   void CoupledLinearNS_sparse::load_session_parameters()
    {
    
    	cout << " Loading ROM parameters ..." << endl;
    
    	ROM_started = 0;
    	ongoing_snapshot_computation = 0;
   
    	load_snapshot_data_from_files = m_session->GetParameter("load_snapshot_data_from_files");
		number_of_snapshots = m_session->GetParameter("number_of_snapshots");
		POD_tolerance = m_session->GetParameter("POD_tolerance");
		if (m_session->DefinesParameter("parameter_space_dimension")) 
		{
			parameter_space_dimension = m_session->GetParameter("parameter_space_dimension");	
		}
		else
		{
			parameter_space_dimension = 1;
		}
		if (m_session->DefinesParameter("debug_mode")) 
		{
			debug_mode = m_session->GetParameter("debug_mode");	
		}
		else
		{
			debug_mode = 1;
		}	
		
		if (m_session->DefinesParameter("do_trafo_check")) 
		{
			do_trafo_check = m_session->GetParameter("do_trafo_check");	
		}
		else
		{
			do_trafo_check = 1;
		}
		if (m_session->DefinesParameter("do_trafo_check_relative_error")) 
		{
			do_trafo_check_relative_error = m_session->GetParameter("do_trafo_check_relative_error");	
		}
		else
		{
			do_trafo_check_relative_error = 1e-9;
		}
		if (m_session->DefinesParameter("compute_smaller_model_errs")) 
		{
			compute_smaller_model_errs = m_session->GetParameter("compute_smaller_model_errs");	
		}
		else
		{
			compute_smaller_model_errs = 0;
		}
		if (m_session->DefinesParameter("use_fine_grid_VV")) 
		{
			use_fine_grid_VV = m_session->GetParameter("use_fine_grid_VV");	
		}
		else
		{
			use_fine_grid_VV = 0;
		}
		if (m_session->DefinesParameter("use_fine_grid_VV_and_load_ref")) 
		{
			use_fine_grid_VV_and_load_ref = m_session->GetParameter("use_fine_grid_VV_and_load_ref");	
		}
		else
		{
			use_fine_grid_VV_and_load_ref = 1;
		}	
		if (m_session->DefinesParameter("use_fine_grid_VV_random")) 
		{
			use_fine_grid_VV_random = m_session->GetParameter("use_fine_grid_VV_random");	
		}
		else
		{
			use_fine_grid_VV_random = 0;
		}
		if (m_session->DefinesParameter("use_sparse_poly")) 
		{
			use_sparse_poly = m_session->GetParameter("use_sparse_poly");
		}
		else
		{
			use_sparse_poly = 0;
		} 
		if (m_session->DefinesParameter("max_sparse_poly_approx_dimension")) 
		{
			max_sparse_poly_approx_dimension = m_session->GetParameter("max_sparse_poly_approx_dimension");
		}
		else
		{
			max_sparse_poly_approx_dimension = 1;
		} 
		if (m_session->DefinesParameter("Geo_trafo_load")) 
		{
			Geo_trafo_load = m_session->GetParameter("Geo_trafo_load");
		}
		else
		{
			Geo_trafo_load = 1;
		} 	
		if (m_session->DefinesParameter("replace_snapshot_with_transformed")) 
		{
			replace_snapshot_with_transformed = m_session->GetParameter("replace_snapshot_with_transformed");
		}
		else
		{
			replace_snapshot_with_transformed = 1;
		} 
	
	
	
		parameter_types = Array<OneD, int> (parameter_space_dimension); 
		parameter_types[0] = m_session->GetParameter("type_para1");
		Nmax = number_of_snapshots;
	   	if (parameter_space_dimension == 1)
		{
			param_vector = Array<OneD, NekDouble> (Nmax);
			for(int i = 0; i < number_of_snapshots; ++i)
			{
				// generate the correct string
				std::stringstream sstm;
				sstm << "param" << i;
				std::string result = sstm.str();
				// const char* rr = result.c_str();
		        param_vector[i] = m_session->GetParameter(result);
	        }
		}
		if (parameter_space_dimension == 2)
		{
			parameter_types[1] = m_session->GetParameter("type_para2");
			general_param_vector = Array<OneD, Array<OneD, NekDouble> > (Nmax);
			number_of_snapshots_dir0 = m_session->GetParameter("number_of_snapshots_dir0");
			number_of_snapshots_dir1 = m_session->GetParameter("number_of_snapshots_dir1");
			int i_all = 0;
			Array<OneD, NekDouble> parameter_point(parameter_space_dimension, 0.0);
			for(int i = 0; i < Nmax; ++i)
			{
				parameter_point = Array<OneD, NekDouble> (parameter_space_dimension, 0.0);
				general_param_vector[i] = parameter_point;
			}
			for(int i0 = 0; i0 < number_of_snapshots_dir0; ++i0)
			{
				// generate the correct string
				std::stringstream sstm;
				sstm << "param" << i0 << "_dir0";
				std::string result = sstm.str();
		        if (i0 == 0)
					start_param_dir0 = m_session->GetParameter(result);
		        if (i0 == number_of_snapshots_dir0-1)
					end_param_dir0 = m_session->GetParameter(result);
				param_vector_dir0[i0] = m_session->GetParameter(result);
				for(int i1 = 0; i1 < number_of_snapshots_dir1; ++i1)
				{
					// generate the correct string
					std::stringstream sstm1;
					sstm1 << "param" << i1 << "_dir1";
					std::string result1 = sstm1.str();
				    if (i1 == 0)
						start_param_dir1 = m_session->GetParameter(result1);
				    if (i1 == number_of_snapshots_dir1-1)
						end_param_dir1 = m_session->GetParameter(result1);
					general_param_vector[i_all][0] = m_session->GetParameter(result);
				    general_param_vector[i_all][1] = m_session->GetParameter(result1);
					i_all = i_all + 1;
	//				general_param_vector[i_all] = Array<OneD, NekDouble> (parameter_space_dimension);
					param_vector_dir1[i1] = m_session->GetParameter(result1);
				}
			}
			for(int i0 = 0; i0 < number_of_snapshots_dir0; ++i0)
			{
				cout << "param_vector_dir0[i0] " << param_vector_dir0[i0] << endl;
			}
			for(int i0 = 0; i0 < number_of_snapshots_dir1; ++i0)
			{
				cout << "param_vector_dir1[i0] " << param_vector_dir1[i0] << endl;
			}
		
		
		}



			for(int i0 = 0; i0 < number_of_snapshots_dir0; ++i0)
			{
				cout << "param_vector_dir0[i0] " << param_vector_dir0[i0] << endl;
			}
			for(int i0 = 0; i0 < number_of_snapshots_dir1; ++i0)
			{
				cout << "param_vector_dir1[i0] " << param_vector_dir1[i0] << endl;
			}




		if (use_fine_grid_VV && (parameter_space_dimension == 1))
		{
			fine_grid_dir0 = m_session->GetParameter("fine_grid_dir0");
			Nmax_VV = fine_grid_dir0;
			fine_general_param_vector = Array<OneD, Array<OneD, NekDouble> >(fine_grid_dir0);
			for(int i = 0; i < fine_grid_dir0; ++i)
			{
				Array<OneD, NekDouble> parameter_point = Array<OneD, NekDouble> (parameter_space_dimension, 0.0);
				fine_general_param_vector[i] = parameter_point;
			}

			// or alternatively overwrite with given random vector
			if (use_fine_grid_VV_random)
			{
				int current_index = 0;
				for (int i1 = 0; i1 < fine_grid_dir0; i1++)
				{
					// generate the correct string
					std::stringstream sstm;
					sstm << "VV_param" << i1 << "_dir0";
					std::string result = sstm.str();
					double param_dir0 = m_session->GetParameter(result);
	

					fine_general_param_vector[current_index][0] = param_dir0;
					current_index++;
 				}
				
			}

		}



/*		if (m_session->DefinesParameter("number_elem_trafo"))
		{
			number_elem_trafo = m_session->GetParameter("number_elem_trafo");
		}
		else
		{
			number_elem_trafo = 5;
		}

		switch(Geo_trafo_load) 
		{
			case 1:
				elements_trafo = set_elem_trafo( number_elem_trafo);
				break;
			case 2:
				elements_trafo = set_elem_trafo_no2( number_elem_trafo);
				break;
			case 3:
				elements_trafo = set_elem_trafo_no3( number_elem_trafo);
				break;
			case 4:
				elements_trafo = set_elem_trafo_no4( number_elem_trafo);
				break;
			case 5:
				elements_trafo = set_elem_trafo_no5( number_elem_trafo);
				break;
			default:
				cout << "error: geo trafo number not known " << endl;				
		}				

*/
	
    	cout << " ... finished loading ROM parameters" << endl;
	

		

    }


    void CoupledLinearNS_sparse::v_InitObject()
    {
		cout << "Hello from CoupledLinearNS_sparse::v_InitObject(void)" << endl;

//        IncNavierStokes::v_InitObject();
        CoupledLinearNS::v_InitObject();


	}


    void CoupledLinearNS_sparse::v_DoSolve(void)
    {
		cout << "Hello from CoupledLinearNS_sparse::v_DoSolve(void)" << endl;

		if (parameter_space_dimension == 1)
		{
        	compute_sparse_poly_approx();
        }
        if (parameter_space_dimension == 2)
		{
        	compute_sparse_poly_approx_2D();
        }

	}

    void CoupledLinearNS_sparse::compute_sparse_poly_approx()
    {



	int sparse_poly_approx_dimension = max_sparse_poly_approx_dimension;
	// L2 error works on the snapshot_x_collection and snapshot_y_collection
	// start sweeping 
	for (int iter_index = 0; iter_index < Nmax; ++iter_index)
	{
		int current_index = iter_index;
		double current_nu;
		double current_geo;
		if (parameter_space_dimension == 1)
		{
			if (parameter_types[0] == 0)
			{
				current_nu = param_vector[current_index];
			}
			if (parameter_types[0] == 1)
			{
				current_geo = param_vector[current_index];
				current_nu = m_kinvis;
			}
		}		
		Array<OneD, NekDouble> interpolant_x(snapshot_x_collection[0].size());
		Array<OneD, NekDouble> interpolant_y(snapshot_y_collection[0].size());
		for (int i = 0; i < snapshot_x_collection[0].size(); ++i)
		{
			interpolant_x[i] = 0.0;
			interpolant_y[i] = 0.0;
		}
		for (int index_interpol_op = 0; index_interpol_op < sparse_poly_approx_dimension; ++index_interpol_op)
		{
			double lagrange_value = lagrange_interp(current_geo, index_interpol_op, sparse_poly_approx_dimension);
			for (int i = 0; i < snapshot_x_collection[0].size(); ++i)	
			{
				interpolant_x[i] += snapshot_x_collection[index_interpol_op][i] * lagrange_value;
				interpolant_y[i] += snapshot_y_collection[index_interpol_op][i] * lagrange_value;
			}
		}
		double rel_L2error = L2norm_abs_error_ITHACA(interpolant_x, interpolant_y, snapshot_x_collection[iter_index], snapshot_y_collection[iter_index]) / L2norm_ITHACA(snapshot_x_collection[iter_index], snapshot_y_collection[iter_index]);
		cout << "rel_L2error at parameter " << current_geo << " is " << rel_L2error << endl;
	}
	Array<OneD, NekDouble> collect_max(sparse_poly_approx_dimension);
	Array<OneD, NekDouble> collect_mean(sparse_poly_approx_dimension);
	double max, mean;
	for (int approx_dim = 1; approx_dim <= sparse_poly_approx_dimension; ++approx_dim)
	{
//		sparse_approx_VV(approx_dim, max, mean);
		cout << "max at dim " << approx_dim << " is " << max << endl;
		cout << "mean at dim " << approx_dim << " is " << mean << endl;
		collect_max[approx_dim-1] = max;
		collect_mean[approx_dim-1] = mean;
	}
	{
	std::stringstream sstm;
	sstm << "sparse_conv_mean.txt";
	std::string sparse_conv_mean = sstm.str();
	const char* outname = sparse_conv_mean.c_str();
	ofstream myfile (outname);
	if (myfile.is_open())
	{
		for (int i0 = 0; i0 < sparse_poly_approx_dimension; i0++)
		{
			myfile << std::setprecision(17) << collect_mean[i0] << "\t";
		}
		myfile.close();
	}
	else cout << "Unable to open file"; 
	}
	{
	std::stringstream sstm;
	sstm << "sparse_conv_max.txt";
	std::string sparse_conv_max = sstm.str();
	const char* outname = sparse_conv_max.c_str();
	ofstream myfile (outname);
	if (myfile.is_open())
	{
		for (int i0 = 0; i0 < sparse_poly_approx_dimension; i0++)
		{
			myfile << std::setprecision(17) << collect_max[i0] << "\t";
		}
		myfile.close();
	}
	else cout << "Unable to open file"; 
	}
    }

	void CoupledLinearNS_sparse::compute_sparse_poly_approx_2D_lower_dim(int approx_dim, double &max_rel_L2, double &mean_rel_L2)
	{
		int sparse_poly_approx_dimension = approx_dim;
		Array<OneD,  Array<OneD, int> > index_set(max_sparse_poly_approx_dimension);
		for (int i=0; i < max_sparse_poly_approx_dimension; ++i)
		{
			index_set[i] = Array<OneD, int>(parameter_space_dimension);
		}
		int added_indices = -1;
		int curr_sum = 0;
		while(added_indices < max_sparse_poly_approx_dimension)
		{
			for(int i=0; i <= curr_sum; ++i)
			{
				++added_indices;
				if (added_indices < max_sparse_poly_approx_dimension)
				{
					int index1 = i;
					int index2 = curr_sum - i;
					index_set[added_indices][0] = index1;
					index_set[added_indices][1] = index2;
				}
			}
			curr_sum++;
		}
		// compute the coefficients
		Array<OneD, Array<OneD, NekDouble> > sparse_poly_coefficients_x(max_sparse_poly_approx_dimension);
		Array<OneD, Array<OneD, NekDouble> > sparse_poly_coefficients_y(max_sparse_poly_approx_dimension);
		for (int i=0; i < max_sparse_poly_approx_dimension; ++i)
		{
			sparse_poly_coefficients_x[i] = Array<OneD, NekDouble> (snapshot_x_collection[0].size());
			sparse_poly_coefficients_y[i] = Array<OneD, NekDouble> (snapshot_x_collection[0].size());
			for(int k=0; k<snapshot_x_collection[0].size(); ++k)
			{
				sparse_poly_coefficients_x[i][k] = 0;
				sparse_poly_coefficients_y[i][k] = 0;
			}
		}
		for (int i=0; i < max_sparse_poly_approx_dimension; ++i)
		{
			// identify the correct snapshot index based on the index_set
			int index_all = index_set[i][1] + number_of_snapshots_dir1 * index_set[i][0]; 
			for(int k=0; k<snapshot_x_collection[0].size(); ++k)
			{
				sparse_poly_coefficients_x[i][k] =  snapshot_x_collection[index_all][k];
				sparse_poly_coefficients_y[i][k] =  snapshot_y_collection[index_all][k];
			}
			for (int j=0; j < i; ++j)
			{
				double lith = lagrange_interp_tensorised_hierarchical(general_param_vector[index_all], index_set[j]);
				for(int k=0; k<snapshot_x_collection[0].size(); ++k)
				{
					sparse_poly_coefficients_x[i][k] -= sparse_poly_coefficients_x[j][k] * lith;
					sparse_poly_coefficients_y[i][k] -= sparse_poly_coefficients_y[j][k] * lith;
				}
			}
		}	
		Array<OneD, NekDouble> collect_L2(Nmax);
		// start sweeping 
		for (int iter_index = 0; iter_index < Nmax; ++iter_index)
		{
			int current_index = iter_index;
			double current_nu;
			double current_geo;
			Array<OneD, NekDouble> current_param(parameter_space_dimension, 0.0);
			current_param = general_param_vector[current_index];
			if (parameter_types[0] == 0)
			{
				current_nu = current_param[0];
				current_geo = current_param[1];
			}
			if (parameter_types[0] == 1)
			{
				current_nu = current_param[1];
				current_geo = current_param[0];
			}
			Array<OneD, NekDouble> interpolant_x(snapshot_x_collection[0].size());
			Array<OneD, NekDouble> interpolant_y(snapshot_y_collection[0].size());
			for (int i = 0; i < snapshot_x_collection[0].size(); ++i)
			{
				interpolant_x[i] = 0.0;
				interpolant_y[i] = 0.0;
			}
			for (int index_interpol_op = 0; index_interpol_op < sparse_poly_approx_dimension; ++index_interpol_op)
			{
				double lith = lagrange_interp_tensorised_hierarchical(general_param_vector[current_index], index_set[index_interpol_op]);
				for (int i = 0; i < snapshot_x_collection[0].size(); ++i)	
				{
					interpolant_x[i] += sparse_poly_coefficients_x[index_interpol_op][i] * lith;
					interpolant_y[i] += sparse_poly_coefficients_y[index_interpol_op][i] * lith;
				}
			}
			double rel_L2error = L2norm_abs_error_ITHACA(interpolant_x, interpolant_y, snapshot_x_collection[iter_index], snapshot_y_collection[iter_index]) / L2norm_ITHACA(snapshot_x_collection[iter_index], snapshot_y_collection[iter_index]);
			collect_L2[iter_index] = rel_L2error;
		}
		mean_rel_L2 = 0;
		max_rel_L2 = 0;
		for (int iter_index = 0; iter_index < Nmax; ++iter_index)
		{
			if (collect_L2[iter_index] > max_rel_L2)
				max_rel_L2 = collect_L2[iter_index];
			mean_rel_L2 += collect_L2[iter_index] / Nmax;
		}
	}

    double CoupledLinearNS_sparse::lagrange_interp(double curr_param, int curr_index, int sparse_poly_approx_dimension)
    {
	double lagrange_value = 1;
	for (int i = 0; i < sparse_poly_approx_dimension; ++i)
	{
		if (parameter_space_dimension == 1)
		{
			if (param_vector[curr_index] != param_vector[i])
				lagrange_value *= (curr_param - param_vector[i]) / (param_vector[curr_index] - param_vector[i]);
		}
	}
	return lagrange_value;
    }



    double CoupledLinearNS_sparse::L2norm_abs_error_ITHACA( Array< OneD, NekDouble > component1_x, Array< OneD, NekDouble > component1_y, Array< OneD, NekDouble > component2_x, Array< OneD, NekDouble > component2_y )
    {
	Array< OneD, NekDouble > x_difference(component1_x.size());
	Array< OneD, NekDouble > y_difference(component1_y.size());
	for (int i = 0; i < component1_y.size(); ++i)
	{
		x_difference[i] = component1_x[i] - component2_x[i];
		y_difference[i] = component1_y[i] - component2_y[i];
	}
	double result = L2norm_ITHACA(x_difference, y_difference);
	return result;
    }

    double CoupledLinearNS_sparse::Linfnorm_abs_error_ITHACA( Array< OneD, NekDouble > component1_x, Array< OneD, NekDouble > component1_y, Array< OneD, NekDouble > component2_x, Array< OneD, NekDouble > component2_y )
    {
	Array< OneD, NekDouble > x_difference(component1_x.size());
	Array< OneD, NekDouble > y_difference(component1_y.size());
	for (int i = 0; i < component1_y.size(); ++i)
	{
		x_difference[i] = component1_x[i] - component2_x[i];
		y_difference[i] = component1_y[i] - component2_y[i];
	}
	double result = Linfnorm_ITHACA(x_difference, y_difference);
	return result;
    }    

    double CoupledLinearNS_sparse::L2norm_ITHACA( Array< OneD, NekDouble > component_x, Array< OneD, NekDouble > component_y )
    {
	// the input comes in phys coords
        NekDouble L2norm = -1.0;
//	cout << "m_NumQuadPointsError " << m_NumQuadPointsError << endl; // should be 0
//	cout << "component_x.size() " << component_x.size() << endl; // should be nphys
        if (m_NumQuadPointsError == 0)
        {
		double L2norm_x = m_fields[0]->L2(component_x);
		double L2norm_y = m_fields[1]->L2(component_y);
		L2norm = sqrt( L2norm_x*L2norm_x + L2norm_y*L2norm_y );
        }
	return L2norm;
    }

    double CoupledLinearNS_sparse::Linfnorm_ITHACA( Array< OneD, NekDouble > component_x, Array< OneD, NekDouble > component_y )
    {
	// the input comes in phys coords
        NekDouble Linfnorm = -1.0;
//	cout << "m_NumQuadPointsError " << m_NumQuadPointsError << endl; // should be 0
//	cout << "component_x.size() " << component_x.size() << endl; // should be nphys
        if (m_NumQuadPointsError == 0)
        {
		double Linfnorm_x = m_fields[0]->L2(component_x);
		double Linfnorm_y = m_fields[1]->L2(component_y);
		Linfnorm = max(Linfnorm_x, Linfnorm_y);
        }
	return Linfnorm;
    }



    
    double CoupledLinearNS_sparse::lagrange_interp_tensorised_hierarchical(Array<OneD, double> curr_param, Array<OneD, int> curr_index)
    {
	double lagrange_value = 1;




//	cout << "curr_param[0] " << curr_param[0] << " curr_param[1] " << curr_param[1] << endl;
//	cout << "curr_index[0] " << curr_index[0] << " curr_index[1] " << curr_index[1] << endl; 
	for (int i = 0; i < curr_index[0]; ++i)
	{
		lagrange_value *= (curr_param[0] - param_vector_dir0[i]) / (param_vector_dir0[curr_index[0]] - param_vector_dir0[i]);
//		cout << "(curr_param[0] - param_vector_dir0[i]) / (param_vector_dir0[curr_index[0]] - param_vector_dir0[i]); " << (curr_param[0] - param_vector_dir0[i]) / (param_vector_dir0[curr_index[0]] - param_vector_dir0[i]) << endl;
//		cout << "param_vector_dir0[curr_index[0]] " << param_vector_dir0[curr_index[0]] << endl;
//		cout << "curr_index[0] " << curr_index[0]	 << endl;
//		cout << "param_vector_dir0[i] " << param_vector_dir0[i] << endl;
	}
	if (parameter_space_dimension > 1)
	{
		for (int i = 0; i < curr_index[1]; ++i)
		{
			lagrange_value *= (curr_param[1] - param_vector_dir1[i]) / (param_vector_dir1[curr_index[1]] - param_vector_dir1[i]);
//			cout << "(curr_param[1] - param_vector_dir1[i]) / (param_vector_dir1[curr_index[1]] - param_vector_dir1[i]); " << (curr_param[1] - param_vector_dir1[i]) / (param_vector_dir1[curr_index[1]] - param_vector_dir1[i]) << endl;
		}
	}
//	cout << "curr_param[0] " << curr_param[0] << endl;
//	cout << "curr_param[1] " << curr_param[1] << endl;  
	return lagrange_value;
    }    



    void CoupledLinearNS_sparse::compute_sparse_poly_approx_2D()
    {


			for(int i0 = 0; i0 < number_of_snapshots_dir0; ++i0)
			{
				cout << "param_vector_dir0[i0] " << param_vector_dir0[i0] << endl;
			}
			for(int i0 = 0; i0 < number_of_snapshots_dir1; ++i0)
			{
				cout << "param_vector_dir1[i0] " << param_vector_dir1[i0] << endl;
			}

	int sparse_poly_approx_dimension = max_sparse_poly_approx_dimension;
	// L2 error works on the snapshot_x_collection and snapshot_y_collection
	
	// need to decide the grid rule 
	Array<OneD,  Array<OneD, int> > index_set(max_sparse_poly_approx_dimension);
	
	for (int i=0; i < max_sparse_poly_approx_dimension; ++i)
	{
		index_set[i] = Array<OneD, int>(parameter_space_dimension);
	}
	// 0,0
	// 1,0
	// 0,1
	// 1,1
	// 2,0
	// 0,2
	// 3,0 ..
	int added_indices = -1;
	int curr_sum = 0;
	while(added_indices < max_sparse_poly_approx_dimension)
	{
		for(int i=0; i <= curr_sum; ++i)
		{
			++added_indices;
			if (added_indices < max_sparse_poly_approx_dimension)
			{
				int index1 = i;
				int index2 = curr_sum - i;
//				cout << "index1 " << index1 << " index2 " << index2 << endl; 				
				index_set[added_indices][0] = index1;
				index_set[added_indices][1] = index2;
			}
		
		}
		curr_sum++;
	}
	
	
	
	// compute the coefficients
	Array<OneD, Array<OneD, NekDouble> > sparse_poly_coefficients_x(max_sparse_poly_approx_dimension);
	Array<OneD, Array<OneD, NekDouble> > sparse_poly_coefficients_y(max_sparse_poly_approx_dimension);
	for (int i=0; i < max_sparse_poly_approx_dimension; ++i)
	{
		sparse_poly_coefficients_x[i] = Array<OneD, NekDouble> (snapshot_x_collection[0].size());
		sparse_poly_coefficients_y[i] = Array<OneD, NekDouble> (snapshot_x_collection[0].size());
		for(int k=0; k<snapshot_x_collection[0].size(); ++k)
		{
			sparse_poly_coefficients_x[i][k] = 0;
			sparse_poly_coefficients_y[i][k] = 0;
		}
	}
	for (int i=0; i < max_sparse_poly_approx_dimension; ++i)
	{
		// identify the correct snapshot index based on the index_set
		int index_all = index_set[i][1] + number_of_snapshots_dir1 * index_set[i][0]; 
		for(int k=0; k<snapshot_x_collection[0].size(); ++k)
		{
			sparse_poly_coefficients_x[i][k] =  snapshot_x_collection[index_all][k];
			sparse_poly_coefficients_y[i][k] =  snapshot_y_collection[index_all][k];
		}
		
		for (int j=0; j < i; ++j)
		{
//			cout << "current j: " << j << " current i: " << i << endl;
			double lith = lagrange_interp_tensorised_hierarchical(general_param_vector[index_all], index_set[j]);
//			cout << "lith " << lith << endl;
			for(int k=0; k<snapshot_x_collection[0].size(); ++k)
			{
				sparse_poly_coefficients_x[i][k] -= sparse_poly_coefficients_x[j][k] * lith;
				sparse_poly_coefficients_y[i][k] -= sparse_poly_coefficients_y[j][k] * lith;
			}
		}
		



	}	

	Array<OneD, NekDouble> collect_L2(Nmax);
	
	// start sweeping 
	for (int iter_index = 0; iter_index < Nmax; ++iter_index)
	{
		int current_index = iter_index;
		double current_nu;
		double current_geo;
		Array<OneD, NekDouble> current_param(parameter_space_dimension, 0.0);
		current_param = general_param_vector[current_index];
		if (parameter_types[0] == 0)
		{
			current_nu = current_param[0];
			current_geo = current_param[1];
		}
		if (parameter_types[0] == 1)
		{
			current_nu = current_param[1];
			current_geo = current_param[0];
		}
		Array<OneD, NekDouble> interpolant_x(snapshot_x_collection[0].size());
		Array<OneD, NekDouble> interpolant_y(snapshot_y_collection[0].size());
		for (int i = 0; i < snapshot_x_collection[0].size(); ++i)
		{
			interpolant_x[i] = 0.0;
			interpolant_y[i] = 0.0;
		}
		for (int index_interpol_op = 0; index_interpol_op < sparse_poly_approx_dimension; ++index_interpol_op)
		{
			double lith = lagrange_interp_tensorised_hierarchical(general_param_vector[current_index], index_set[index_interpol_op]);
//			cout << "lith sweep " << lith << endl;
			for (int i = 0; i < snapshot_x_collection[0].size(); ++i)	
			{
				interpolant_x[i] += sparse_poly_coefficients_x[index_interpol_op][i] * lith;
				interpolant_y[i] += sparse_poly_coefficients_y[index_interpol_op][i] * lith;
//				interpolant_x[i] += sparse_poly_coefficients_x[index_interpol_op][i];
//				interpolant_y[i] += sparse_poly_coefficients_y[index_interpol_op][i];
			}
		}
		cout << " interpolant_x[0] " << interpolant_x[0] << endl;
		cout << " snapshot_x_collection[0][0] " << snapshot_x_collection[0][0] << endl;
		double rel_L2error = L2norm_abs_error_ITHACA(interpolant_x, interpolant_y, snapshot_x_collection[iter_index], snapshot_y_collection[iter_index]) / L2norm_ITHACA(snapshot_x_collection[iter_index], snapshot_y_collection[iter_index]);
		cout << "rel_L2error at 2D parameter geo " << current_geo << " nu " << current_nu << " is " << rel_L2error << endl;
		collect_L2[iter_index] = rel_L2error;
	}
	double mean_rel_L2 = 0;
	double max_rel_L2 = 0;
	for (int iter_index = 0; iter_index < Nmax; ++iter_index)
	{
		if (collect_L2[iter_index] > max_rel_L2)
			max_rel_L2 = collect_L2[iter_index];
		mean_rel_L2 += collect_L2[iter_index] / Nmax;
	}
	cout << "mean_rel_L2 " << mean_rel_L2 << " max_rel_L2 " << max_rel_L2 << endl;
	Array<OneD, NekDouble> collect_max(sparse_poly_approx_dimension);
	Array<OneD, NekDouble> collect_mean(sparse_poly_approx_dimension);
	double max, mean;
	for (int approx_dim = 1; approx_dim <= sparse_poly_approx_dimension; ++approx_dim)
	{
		compute_sparse_poly_approx_2D_lower_dim(approx_dim, max_rel_L2, mean_rel_L2);
		cout << "max at dim " << approx_dim << " is " << max_rel_L2 << endl;
		cout << "mean at dim " << approx_dim << " is " << mean_rel_L2 << endl;
		collect_max[approx_dim-1] = max_rel_L2;
		collect_mean[approx_dim-1] = mean_rel_L2;
	}
	{
	std::stringstream sstm;
	sstm << "sparse_conv_mean.txt";
	std::string sparse_conv_mean = sstm.str();
	const char* outname = sparse_conv_mean.c_str();
	ofstream myfile (outname);
	if (myfile.is_open())
	{
		for (int i0 = 0; i0 < sparse_poly_approx_dimension; i0++)
		{
			myfile << std::setprecision(17) << collect_mean[i0] << "\t";
		}
		myfile.close();
	}
	else cout << "Unable to open file"; 
	}
	{
	std::stringstream sstm;
	sstm << "sparse_conv_max.txt";
	std::string sparse_conv_max = sstm.str();
	const char* outname = sparse_conv_max.c_str();
	ofstream myfile (outname);
	if (myfile.is_open())
	{
		for (int i0 = 0; i0 < sparse_poly_approx_dimension; i0++)
		{
			myfile << std::setprecision(17) << collect_max[i0] << "\t";
		}
		myfile.close();
	}
	else cout << "Unable to open file"; 
	}
    }



    void CoupledLinearNS_sparse::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
		cout << "Hello from CoupledLinearNS_sparse::v_GenerateSummary" << endl;
	}

        
       

    bool CoupledLinearNS_sparse::v_NegatedOp(void)
    {
		cout << "Hello from CoupledLinearNS_sparse::v_NegatedOp" << endl;
	}


    void CoupledLinearNS_sparse::v_TransCoeffToPhys()
    {
		cout << "Hello from CoupledLinearNS_sparse::v_TransCoeffToPhys" << endl;
	}

        

    void CoupledLinearNS_sparse::v_TransPhysToCoeff()
    {
		cout << "Hello from CoupledLinearNS_sparse::v_TransPhysToCoeff" << endl;
	}

       



    void CoupledLinearNS_sparse::v_Output()
    {
		cout << "Hello from CoupledLinearNS_sparse::v_Output" << endl;
	}

        

    int CoupledLinearNS_sparse::v_GetForceDimension()
    {
		cout << "Hello from CoupledLinearNS_sparse::v_GetForceDimension" << endl;
	}




}
