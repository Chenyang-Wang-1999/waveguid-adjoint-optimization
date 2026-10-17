#include "EM_mesh_flagging.h"
#include <libmesh/mesh.h>
#include <libmesh/elem.h>
#include <libmesh/boundary_info.h>

using namespace libMesh;

EMMeshFlagging::EMMeshFlagging(MeshBase * mesh,
                        std::vector<subdomain_id_type> &r_id,
                        std::vector<double> &r_index,
                        double default_index,
                        std::vector<subdomain_id_type> &p_id,
                        double port_size) : regions_id(r_id),regions_index(r_index),port_id(p_id)
{
    this->mesh = mesh;
    //this->regions_id = regions_id;
    //this->regions_index = regions_index;
    this->default_index = default_index;
    this->N_regions = regions_id.size();
    //this->port_id = port_id;
    this->port_size = port_size;
    this->N_ports = port_id.size();
}

void EMMeshFlagging::set_freq(double freq)
{
    this->freq = freq;
}

void EMMeshFlagging::flag_elements()
{
    // refine the mesh accroding to the frequency and the permitivity
    // set parameters
    // std::cout << regions_index[0]<<'\n';
    unsigned int dim = mesh->mesh_dimension();
    BoundaryInfo mesh_bd = mesh->get_boundary_info();
    // traverse all elements
    for(Elem* elem : mesh->active_local_element_ptr_range())
    {
        double max_elem_vol;
        double index;
        bool is_default = true;
        bool is_port = false;
        // which domain
        // std::cout << "N regions:"<<N_regions<<'\n';
        for(unsigned int curr_domain = 0; curr_domain < N_regions; curr_domain ++)
        {
            subdomain_id_type curr_id = regions_id[curr_domain];
            /*
            std::cout << "elem subdomain id:" << elem->subdomain_id() <<'\n';
            std::cout << "current id:" << curr_id << '\n';
            //*/
            if(elem->subdomain_id() == curr_id)
            {

                is_default = false;
                index = regions_index[curr_domain];
                /*
                std::cout << "This subdomain is "<<curr_id <<'\n';
                std::cout << "Index is:" << index <<'\n';
                //*/
                break;
            }
        }
        if(is_default)
        {
            index = default_index;
        }

        // port refinement
        for(unsigned int curr_port = 0; curr_port < N_ports; curr_port ++)
        {
            subdomain_id_type curr_id = port_id[curr_port];
            for(auto side: elem->side_index_range())
            {
                if(mesh_bd.has_boundary_id(elem, side, curr_id))
                {
                    is_port = true;
                    std::unique_ptr<Elem> side_elem = elem->side_ptr(side);
                    // std::cout << "current port size" << side_elem -> volume() << '\n';
                    if(side_elem->volume()>pow(port_size,dim-1))
                    {
                        elem->set_refinement_flag(Elem::REFINE);
                    }
                    else
                    {
                        elem->set_refinement_flag(Elem::DO_NOTHING);
                    }
                    break;
                }
            }
            if(is_port)
            {
                break;
            }

        }


        // get wavelength
        double curr_wl = 1/freq/index;
        max_elem_vol = pow(curr_wl/16, dim);
        /* debug
        std::cout << "Current index: "<< index <<'\n';
        std::cout << "Element volume: " << elem->volume() <<'\n';
        std::cout << "Max volume: " << max_elem_vol <<'\n';
        //*/
        // set refinement
        if(elem->volume() > max_elem_vol)
        {
            elem->set_refinement_flag(Elem::REFINE);
        }
        else
        {
            if(!is_port)
            {
                elem->set_refinement_flag(Elem::DO_NOTHING);
            }
        }
    }
}
