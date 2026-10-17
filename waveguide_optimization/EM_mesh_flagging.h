#include <libmesh/mesh_refinement.h>
#include <libmesh/mesh.h>
#include <vector>

#ifndef EM_MESH_FLAGGING_INCLUDED
#define EM_MESH_FLAGGING_INCLUDED

using namespace libMesh;

class EMMeshFlagging: public MeshRefinement::ElementFlagging
{
    public:
        EMMeshFlagging(MeshBase * mesh,
                        std::vector<subdomain_id_type> &r_id,
                        std::vector<double> &r_index,
                        double default_index,
                        std::vector<subdomain_id_type> &p_id,
                        double port_size);
        void flag_elements();
        void set_freq(double freq);

    private:
        std::vector<subdomain_id_type> &regions_id;
        std::vector<double> &regions_index;
        MeshBase * mesh;
        double default_index;
        unsigned int N_regions;
        double freq;
        std::vector<subdomain_id_type> &port_id;
        double port_size;
        unsigned int N_ports;
};

#endif // EM_MESH_FLAGGING_INCLUDED
