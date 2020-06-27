function build_dfd_net_rw_v3(filename, pso_member, net_description)

    con_idx = 1;
    bn_idx = 1;
    act_idx = 1;
    tag_idx = 1;
    concat_idx = 1;
    cond_idx = 1;
    end_cnt = 0;

    input_type = [2,1;3,1;3,2;3,3;6,1;6,2;6,3];
    
    fprintf('// Auto generated c++ header file for PSO testing\n');
    fprintf('// Iteration Number: %d\n', pso_member.iteration);
    fprintf('// Population Number: %d\n\n', pso_member.number);

    % print the net headers
    fprintf('#ifndef NET_DEFINITION_H\n');
    fprintf('#define NET_DEFINITION_H\n\n');
    fprintf('#include <cstdint>\n');
    fprintf('#include "dlib/dnn.h"\n#include "dlib/dnn/core.h"\n\n');
    fprintf('#include "gorgon_capture.h"\n\n');

    fprintf('extern const uint32_t img_depth = %d;\nextern const uint32_t secondary = %d;\n',input_type(pso_member.input,:));
    fprintf('std::vector<std::pair<uint64_t, uint64_t>> crop_sizes = {{%d, %d}, {38, 148}};\n\n',pso_member.crop_size(1)*18,pso_member.crop_size(2)*6);
    

    fprintf('typedef struct{\n    uint32_t iteration;\n    uint32_t pop_num;\n} pso_struct;\n');
    fprintf('pso_struct pso_info = {%d, %d};\n\n', pso_member.iteration, pso_member.number);

    fprintf('template <long num_filters, typename SUBNET> using con2d = dlib::con<num_filters, 2, 2, 2, 2, SUBNET>;\n');
    fprintf('template <long num_filters, typename SUBNET> using con33d33 = dlib::con<num_filters, 3, 3, 3, 3, SUBNET>;\n');
    fprintf('template <long num_filters, typename SUBNET> using con32d32 = dlib::con<num_filters, 3, 2, 3, 2, SUBNET>;\n');
    fprintf('template <long num_filters, typename SUBNET> using con22d21 = dlib::con<num_filters, 2, 1, 2, 1, SUBNET>;\n\n');

    fprintf('template <long num_filters, typename SUBNET> using cont2u = dlib::cont<num_filters, 2, 2, 2, 2, SUBNET>;\n\n');

    fprintf('template <typename SUBNET> using DTO_0 = dlib::add_tag_layer<200, SUBNET>;\n');
    fprintf('template <typename SUBNET> using DTI_0 = dlib::add_tag_layer<201, SUBNET>;\n');
    fprintf('template <typename SUBNET> using DTO_1 = dlib::add_tag_layer<202, SUBNET>;\n');
    fprintf('template <typename SUBNET> using DTI_1 = dlib::add_tag_layer<203, SUBNET>;\n');
    fprintf('template <typename SUBNET> using DTO_2 = dlib::add_tag_layer<204, SUBNET>;\n');
    fprintf('template <typename SUBNET> using DTI_2 = dlib::add_tag_layer<205, SUBNET>;\n\n');

    fprintf('using dfd_net_type = dlib::loss_multiclass_log_per_pixel<\n');

    for idx=1:numel(net_description.net_structure)

        if(strcmp(net_description.net_structure{idx}, 'con'))
            fprintf('    dlib::con<%d, %d, %d, 1, 1, \n', pso_member.con(con_idx,1),2*pso_member.con(con_idx,2)+1,2*pso_member.con(con_idx,3)+1); 
            con_idx = con_idx + 1;
        elseif(strcmp(net_description.net_structure{idx}, 'act'))
            fprintf('    %s', net_description.activations{pso_member.act(act_idx)});
            act_idx = act_idx + 1;
        elseif(strcmp(net_description.net_structure{idx}, 'bn'))
            if(pso_member.bn(bn_idx) == 1)
                fprintf('    dlib::bn_con<');
            else
                end_cnt = end_cnt + 1;
            end
            bn_idx = bn_idx + 1;
        elseif(strcmp(net_description.net_structure{idx}, 'cond'))
            fprintf('    %s%d,', net_description.cond{cond_idx}, pso_member.con(con_idx,1)); 
            cond_idx = cond_idx + 1;
            con_idx = con_idx + 1;
        elseif(strcmp(net_description.net_structure{idx}, 'cont'))
            fprintf('    dlib::cont<%d, %d, %d, 2, 2, \n', pso_member.con(con_idx,1),2*pso_member.con(con_idx,2)+1,2*pso_member.con(con_idx,3)+1) 
            con_idx = con_idx + 1;
        elseif(strcmp(net_description.net_structure{idx}, 'cont2u'))
            fprintf('    dlib::cont<%d, 2, 2, 2, 2, \n', pso_member.con(con_idx,1)) 
            con_idx = con_idx + 1;
        elseif(strcmp(net_description.net_structure{idx}, 'concat'))
            fprintf('    %s\n', net_description.concats{concat_idx}) 
            concat_idx = concat_idx + 1;        
        elseif(strcmp(net_description.net_structure{idx}, 'add_prev1'))
            fprintf('    dlib::add_prev1<');
        elseif(strcmp(net_description.net_structure{idx}, 'tag1'))
            fprintf('    dlib::tag1<');
        elseif(strcmp(net_description.net_structure{idx}, 'tags'))
            fprintf('    %s',net_description.tags{tag_idx});
            tag_idx = tag_idx + 1;
        end    

    end
    
    e = repmat('>',1,idx-end_cnt+1);
    fprintf('    dlib::input<std::array<dlib::matrix<uint16_t>, img_depth>>\n');
    fprintf('    %s;\n\n',e);

    fprintf('//---------------------------------------------------------\n\n');

    fprintf('inline std::ostream& operator<< (std::ostream& out, const pso_struct& item)\n{\n')
    fprintf('    out << item.iteration << ", " << item.pop_num;\n');
    fprintf('    return out;\n}\n\n');

    fprintf('#endif \n');

end