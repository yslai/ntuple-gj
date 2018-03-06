#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <H5Cpp.h>

// This is chosen to be the CPU L2 cache size, which should exceed 512
// kB for many years now
#ifndef HDF5_DEFAULT_CACHE
#define HDF5_DEFAULT_CACHE (512 * 1024)
#endif // HDF5_DEFAULT_CACHE

#ifndef HDF5_USE_DEFLATE
#define HDF5_USE_DEFLATE
#endif // HDF5_USE_DEFLATE

// Rank of the tensor written in HDF5, rank 3 being (index_event,
// index_track, index_properties)
#define RANK 3

void find_ntrack_ncluster_max(char *argv_first[], char *argv_last[], UInt_t &ntrack_max, UInt_t &ncluster_max)
{
 
   for (char **p = argv_first; p != argv_last; p++) {
        // Cautious opening of the TTree, capturing all modes of
        // failure, and keep the TDirectoryFile (to be deleted later)
        // to avoid memory leak
        TFile *file = TFile::Open(*p);

        if (file == NULL) {
            continue;
        }

         TDirectoryFile *df = dynamic_cast<TDirectoryFile *>
             (file->Get("AliAnalysisTaskNTGJ"));

         if (df == NULL) {
	       fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "Cannot open TFile");
             continue;
         }

         TTree *hi_tree = dynamic_cast<TTree *>
             (df->Get("_tree_event"));

	 if (hi_tree == NULL) {
	   fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "Cannot open _tree_event");
	   continue;
        }

        // Get the maximum of the "ntrack"

        UInt_t ntrack;
        UInt_t ncluster;

        hi_tree->SetBranchAddress("ntrack", &ntrack);
        hi_tree->SetBranchAddress("ncluster", &ncluster);

        for (Long64_t i = 0; i < hi_tree->GetEntries(); i++) {
            hi_tree->GetEntry(i);
            ntrack_max = std::max(ntrack_max, ntrack);
	    ncluster_max = std::max(ncluster_max, ncluster);
	    fprintf(stderr, "\r%s:%d: %llu", __FILE__, __LINE__, i);
	}
	fprintf(stderr, "\n");
        // Fully delete everything

        hi_tree->Delete();
	delete df;
        file->Close();
        delete file;
    }

}

void write_track_cluster(H5::DataSet &track_data_set, H5::DataSet &cluster_data_set, 
			 hsize_t *offset, const hsize_t *track_dim_extend, 
			 const hsize_t *cluster_dim_extend, const UInt_t ntrack_max,
			 const UInt_t ncluster_max, char *argv_first[], char *argv_last[])
{
    for (char **p = argv_first; p != argv_last; p++) {
        TFile *file = TFile::Open(*p);

        if (file == NULL) {
            continue;
        }

        TDirectoryFile *df = dynamic_cast<TDirectoryFile *>
            (file->Get("AliAnalysisTaskNTGJ"));

        if (df == NULL) {
            continue;
        }

        TTree *hi_tree = dynamic_cast<TTree *>
            (df->Get("_tree_event"));

        if (hi_tree == NULL) {
            continue;
        }

        UInt_t ntrack;
        std::vector<Float_t> track_e(ntrack_max, NAN);
        std::vector<Float_t> track_pt(ntrack_max, NAN);
        std::vector<Float_t> track_eta(ntrack_max, NAN);
        std::vector<Float_t> track_phi(ntrack_max, NAN);
        std::vector<UChar_t> track_quality(ntrack_max, NAN);
        std::vector<Float_t> track_eta_emcal(ntrack_max, NAN);
        std::vector<Float_t> track_phi_emcal(ntrack_max, NAN);

	UInt_t ncluster;
	std::vector<Float_t> cluster_e(ncluster_max, NAN);
	std::vector<Float_t> cluster_pt(ncluster_max, NAN);
	std::vector<Float_t> cluster_eta(ncluster_max, NAN);
	std::vector<Float_t> cluster_phi(ncluster_max, NAN);
	std::vector<Float_t> cluster_e_cross(ncluster_max, NAN);
	//std::vector<std::vector<Float_t> > cluster_s_nphoton(ncluster_max,std::vector <Float_t> (4, NAN) );

        hi_tree->SetBranchAddress("ntrack", &ntrack);
        hi_tree->SetBranchAddress("track_e", &track_e[0]);
        hi_tree->SetBranchAddress("track_pt", &track_pt[0]);
        hi_tree->SetBranchAddress("track_eta", &track_eta[0]);
        hi_tree->SetBranchAddress("track_phi", &track_phi[0]);
        hi_tree->SetBranchAddress("track_quality", &track_quality[0]);
        hi_tree->SetBranchAddress("track_eta_emcal", &track_eta_emcal[0]);
        hi_tree->SetBranchAddress("track_phi_emcal", &track_phi_emcal[0]);

        hi_tree->SetBranchAddress("ncluster", &ncluster);
        hi_tree->SetBranchAddress("cluster_e", &cluster_e[0]);
        hi_tree->SetBranchAddress("cluster_pt", &cluster_pt[0]);
        hi_tree->SetBranchAddress("cluster_eta", &cluster_eta[0]);
        hi_tree->SetBranchAddress("cluster_phi", &cluster_phi[0]);
        hi_tree->SetBranchAddress("cluster_e_cross", &cluster_e_cross[0]);
	//hi_tree->SetBranchAddress("cluster_s_nphoton", &cluster_s_nphoton[0]);

        for (Long64_t i = 0; i < hi_tree->GetEntries(); i++) {
            hi_tree->GetEntry(i);

            std::vector<float> track_data(ntrack_max * 7, NAN);
            std::vector<float> cluster_data(ncluster_max * 5, NAN);

            for (Long64_t j = 0; j < ntrack; j++) {
                // Note HDF5 is always row-major (C-like)
                track_data[j * 7 + 0] = track_e[j];
                track_data[j * 7 + 1] = track_pt[j];
                track_data[j * 7 + 2] = track_eta[j];
                track_data[j * 7 + 3] = track_phi[j];
                track_data[j * 7 + 4] = track_quality[j];
                track_data[j * 7 + 5] = track_eta_emcal[j];
                track_data[j * 7 + 6] = track_phi_emcal[j];
            }

	    for(Long64_t n = 0; n < ncluster; n++){
	        cluster_data[n * 5 + 0] = cluster_e[n];
	        cluster_data[n * 5 + 1] = cluster_pt[n];
	        cluster_data[n * 5 + 2] = cluster_eta[n];
	        cluster_data[n * 5 + 3] = cluster_phi[n];
	        cluster_data[n * 5 + 4] = cluster_e_cross[n];
	        //cluster_data[j * 7 + 5] = cluster_s_nphoton[j][0];
	        //cluster_data[j * 7 + 6] = cluster_s_nphoton[j][1];
	    }

            if (offset == 0) {
                // Writing the first event. The track_data space is already
                // created with space for one event (see when
                // file.createDataSet() was called)
                track_data_set.write(&track_data[0], H5::PredType::NATIVE_FLOAT);
		cluster_data_set.write(&cluster_data[0], H5::PredType::NATIVE_FLOAT);
            }
            else {
                // The new, extended-by-1 dimension
                const hsize_t track_dim_extended[RANK] = {
                    offset[0] + track_dim_extend[0], track_dim_extend[1],
                    track_dim_extend[2]
                };

		const hsize_t cluster_dim_extended[RANK] = {
		  offset[0] + cluster_dim_extend[0], cluster_dim_extend[1],
		  cluster_dim_extend[2]
                };


                // Extend to the new dimension
                track_data_set.extend(track_dim_extended);
                cluster_data_set.extend(cluster_dim_extended);

                // Select the hyperslab that only encompass the
                // difference from extending the data space (i.e. the
                // new event, but offsetted at the existing event)
                H5::DataSpace track_file_space = track_data_set.getSpace();
		H5::DataSpace cluster_file_space = cluster_data_set.getSpace();

                track_file_space.selectHyperslab(
                    H5S_SELECT_SET, track_dim_extend, offset);

                cluster_file_space.selectHyperslab(
                    H5S_SELECT_SET, cluster_dim_extend, offset);

                // The memory space is the difference only (i.e. also
                // the new event, but at offset 0)
                H5::DataSpace track_memory_space(RANK, track_dim_extend, NULL);
                H5::DataSpace cluster_memory_space(RANK, cluster_dim_extend, NULL);

                // Write the data from memory space to file space
                track_data_set.write(&track_data[0], H5::PredType::NATIVE_FLOAT,
                               track_memory_space, track_file_space);

                cluster_data_set.write(&cluster_data[0], H5::PredType::NATIVE_FLOAT,
                               cluster_memory_space, cluster_file_space);

            }
            offset[0]++;

            if (offset[0] % 1000 == 0) {
                fprintf(stderr, "%s:%d: %llu / %lld\n", __FILE__,
                        __LINE__, offset[0],
                        hi_tree->GetEntries());
            }
        }

        hi_tree->Delete();
        delete df;
        file->Close();
        delete file;
    }
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        exit(EXIT_FAILURE);
    }

    //UInt_t ntrack_max = 1927 ; 
    
    UInt_t ntrack_max = 0;
    UInt_t ncluster_max = 0;

    find_ntrack_ncluster_max(argv + 1, argv + argc - 1, ntrack_max, ncluster_max);
    fprintf(stderr, "%sf:%d: ntrack_max = %u, ncluster_max = %u\n", __FILE__, __LINE__, ntrack_max, ncluster_max);

    // Access mode H5F_ACC_TRUNC truncates any existing file, while
    // not throwing any exception (unlike H5F_ACC_RDWR)
    H5::H5File file(argv[argc - 1], H5F_ACC_TRUNC);

    // How many properties per track is written
    static const size_t track_row_size = 7;
    static const size_t cluster_row_size = 5;
    //easier to just make cluster use same # properties, reuse dim_extend


    // The tensor dimension increment for each new event
    hsize_t track_dim_extend[RANK] = { 1, ntrack_max, track_row_size };
    hsize_t cluster_dim_extend[RANK] = { 1, ncluster_max, cluster_row_size };

    // The maximum tensor dimension, for unlimited number of events
    hsize_t track_dim_max[RANK] = { H5S_UNLIMITED, ntrack_max, track_row_size };
    hsize_t cluster_dim_max[RANK] = { H5S_UNLIMITED, ncluster_max, cluster_row_size };

    // The extensible HDF5 data space
	H5::DataSpace track_data_space(RANK, track_dim_extend, track_dim_max);
	H5::DataSpace cluster_data_space(RANK, cluster_dim_extend, cluster_dim_max);
	//might need two data_spaces: track_data_space & cluster_data_space


    // To enable zlib compression (there will be many NANs) and
    // efficient chunking (splitting of the tensor into contingous
    // hyperslabs), a HDF5 property list is needed
    H5::DSetCreatPropList track_property = H5::DSetCreatPropList();
    H5::DSetCreatPropList cluster_property = H5::DSetCreatPropList();

#ifdef HDF5_USE_DEFLATE
    // Check for zlib (deflate) availability and enable only if
    // present
    if (!H5Zfilter_avail(H5Z_FILTER_DEFLATE)) {
        fprintf(stderr, "%s:%d: warning: deflate filter not "
                "available\n", __FILE__, __LINE__);
    }
    else {
        unsigned int filter_info;

        H5Zget_filter_info(H5Z_FILTER_DEFLATE, &filter_info);
        if (!(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED)) {
            fprintf(stderr, "%s:%d: warning: deflate filter not "
                    "available for encoding\n", __FILE__, __LINE__);
        }
        else {
            track_property.setDeflate(1);
	    cluster_property.setDeflate(1);
        }
    }
#endif // HDF5_USE_DEFLATE

    // Activate chunking, while observing the HDF5_DEFAULT_CACHE being
    // the CPU L2 cache size
    hsize_t track_dim_chunk[RANK] = {
        std::max(static_cast<unsigned long long>(1),
                 HDF5_DEFAULT_CACHE /
                 std::max(static_cast<unsigned long long>(1),
                          track_dim_extend[1] * track_dim_extend[2] *
                          sizeof(float))),
        track_dim_extend[1],
        track_dim_extend[2]
    };

    hsize_t cluster_dim_chunk[RANK] = {
        std::max(static_cast<unsigned long long>(1),
                 HDF5_DEFAULT_CACHE /
                 std::max(static_cast<unsigned long long>(1),
                          cluster_dim_extend[1] * cluster_dim_extend[2] *
                          sizeof(float))),
        cluster_dim_extend[1],
        cluster_dim_extend[2]
    };

    track_property.setChunk(RANK, track_dim_chunk);
    cluster_property.setChunk(RANK, cluster_dim_chunk);

    // Create the data set, which will have space for the first event
    H5::DataSet track_data_set =
        file.createDataSet("track", H5::PredType::NATIVE_FLOAT,
                           track_data_space, track_property);

    H5::DataSet cluster_data_set =
      file.createDataSet("cluster", H5::PredType::NATIVE_FLOAT,
			 cluster_data_space, cluster_property);

    hsize_t offset[RANK] = {0, 0, 0};

    write_track_cluster(track_data_set, cluster_data_set, offset, track_dim_extend, cluster_dim_extend, 
			ntrack_max, ncluster_max, argv + 1, argv + argc - 1);

    file.close();

    return EXIT_SUCCESS;
}
