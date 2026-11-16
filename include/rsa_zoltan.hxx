//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

namespace Zoltan_manager
{
template<typename T>
  struct LBManager
  {
    Zoltan_Struct* zz;
    
    Manager(MPI_COMM& comm)
    {
      zz = Zoltan_Create(comm);
    }

    void set_method(std::string method = "RCB")
    {
			// General parameters 
			//  Zoltan_Set_Param(zz, "LB_METHOD", "RANDOM");    /* Zoltan method: "BLOCK" */
			//  Zoltan_Set_Param(zz, "LB_METHOD", "CYCLIC");    /* Zoltan method: "BLOCK" */
			//  Zoltan_Set_Param(zz, "LB_METHOD", "BLOCK");    /* Zoltan method: "BLOCK" */
			//  Zoltan_Set_Param(zz, "LB_METHOD", "HSFC");    /* Zoltan method: "BLOCK" */
			//  Zoltan_Set_Param(zz, "LB_METHOD", "RIB");    /* Zoltan method: "BLOCK" */
			Zoltan_Set_Param(zz, "LB_METHOD", method);    /* Zoltan method: "BLOCK" */
			Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");  /* global ID is 1 integer */
			Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");  /* local ID is 1 integer */
			Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1");   /* we omit object weights */
			Zoltan_Set_Param(zz, "AUTO_MIGRATE", "TRUE");

      if(method == "RCB")
      {
			  Zoltan_Set_Param(zz, "RCB_RECTILINEAR_BLOCKS","1");
      }
		}

		template<typename T>
			get_new_partition(T& grid)
			{
				// set zoltan query functions depending on the namespace used
				Zoltan_Set_Num_Obj_Fn(zz, get_number_of_objects, &grid);
				Zoltan_Set_Obj_List_Fn(zz, get_object_list, &grid);
				Zoltan_Set_Num_Geom_Fn(zz, fdim, &grid);
				Zoltan_Set_Geom_Multi_Fn(zz, get_geometry_list, &grid);
				Zoltan_Set_Obj_Size_Multi_Fn(zz, user_size_multi_node, &grid);
				Zoltan_Set_Pack_Obj_Multi_Fn(zz, user_pack_multi_node, &grid);
				Zoltan_Set_Unpack_Obj_Multi_Fn(zz, user_unpack_multi_node, &grid);

				// run parition
				int changes, numGidEntries, numLidEntries, numImport, numExport;
				ZOLTAN_ID_PTR importGlobalIds, importLocalIds, exportGlobalIds, exportLocalIds;
				int *importProcs, *importToPart, *exportProcs, *exportToPart;

				int rc = Zoltan_LB_Partition(zz, &changes, &numGidEntries, &numLidEntries,
						&numImport, &importGlobalIds, &importLocalIds, &importProcs, &importToPart,
						&numExport, &exportGlobalIds, &exportLocalIds, &exportProcs, &exportToPart);

				if (rc != ZOLTAN_OK)
				{
					printf("Partitioning failed on process %d\n",rank);
					MPI_Finalize();
					Zoltan_Destroy(&zz);
					exit(0);
				}
			}

		~LBManager() 
		{ 
			Zoltan_Destroy(&zz); 
		}


		static int fdim(void *data, int *ierr){
			*ierr = ZOLTAN_OK;
			return int(3);
		}

		static void get_object_list(void *data, int sizeGID, int sizeLID,
				ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
				int wgt_dim, float *obj_wgts, int *ierr){

			rsa_grid* objs = (graph *)data;
      auto& real = rsa_grid->get_traversal<Traversal::Real>();
			*ierr = ZOLTAN_OK;

			for (int i=0; i< real.size() ; i++)
			{
				globalID[i] = rsa_grid->get_global_id(real[i]);
				localID[i] = real[i];
				obj_wgts[i] = 1; 
			}
			return;
		}

		static int get_number_of_objects(void *data, int *ierr){
			//Catch_Time_Section("zoltan_manager::graph::number_of_objects");
			rsa_grid* objs = (graph *)data;
      auto& real = rsa_grid->get_traversal<Traversal::Real>();
			*ierr = ZOLTAN_OK;

			return real.size();
		}

		static void get_geometry_list(void *data, int sizeGID, int sizeLID,
				int num_obj,
				ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
				int num_dim, double *geom_vec, int *ierr)
		{
			//Catch_Time_Section("zoltan_manager::graph::build_geometry");
			int i;

			rsa_grid* objs = (graph *)data;
      auto& real = rsa_grid->get_traversal<Traversal::Real>();
			*ierr = ZOLTAN_OK;

			assert(num_dim == 3);

			for (i=0;  i < num_obj ; i++)
			{
    
				for(int dim = 0 ; dim < num_dim ; dim++)
				{
					geom_vec[num_dim*i+dim] = (double) (rsa_grid.get_coorvec[i].pos[dim]);
				}
			}

			return;
		}


		static void user_pack_multi_node(void *data,
				int num_gid_entries, int num_lid_entries, int num_ids,
				ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
				int* dest_proc, int* sizes, int* idx,  char *buf, int *ierr)
		{
			Catch_Time_Section("zoltan_manager::graph::pack");
			/* Copy the specified node's data into buffer buf. */
			node *node_buf = (node *) buf;
			graph* objs = (graph *)data;
			auto& vec = objs->get_data();

			for(int i = 0; i < num_ids; i++)
			{
				struct node *node_buf = (struct node *) (buf + idx[i]) ;
				(*node_buf) = vec[local_id[i]];
			}
			/* erase elem */
			std::vector<int> debug(num_ids);
			for(int i = 0; i < num_ids; i++)
			{
				debug[i] = local_id[i];
			}
			std::sort(debug.begin(), debug.end());
			for(int i = debug.size()-1 ; i >= 0 ; i--)
			{
				objs->remove(debug[i]);
			}
			objs->shrink();
			*ierr = ZOLTAN_OK;
		}

		static void user_unpack_multi_node(void *data, int num_gid_entries, int num_ids,
				ZOLTAN_ID_PTR global_id, int* size, int* idx,
				char *buf, int *ierr)
		{
			Catch_Time_Section("zoltan_manager::graph::unpack");
			graph *objs = (graph *)data;
			for(int i = 0; i < num_ids; i++)
			{
				node* node_buf = (node *) (buf+idx[i]);
				objs->add(*node_buf);
			}
			*ierr = ZOLTAN_OK;
		}

	};
};


