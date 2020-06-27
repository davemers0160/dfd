function build_dfd_net_rw(filename, pso_member, net_description)

    con_idx = 1;
    bn_idx = 1;
    act_idx = 1;
    tag_idx = 1;
    concat_idx = 1;
    cond_idx = 1;
    end_cnt = 1;

    input_type = [2,1;3,1;3,2;3,3;6,1;6,2;6,3];
    
    crop_num = [110,110,100,90,65,56,40,28,20,18,16,14,10,8,6,5];
    crop_scale = [4,4];
    %cond_order = [3, 2, 1; 3, 1, 2; 2, 3, 1; 2, 1, 3; 1, 3, 2; 1, 2, 3];
    
    %input_files = {{'dfd_rw_133_128_train_input_v2.txt','dfd_rw_133_128_test_input_v2.txt'};}; 
    input_files = {{'dfd_train_data_sm2.txt','dfd_test_data_sm2.txt'};
                    {'dfd_rw_140_129_train_input_v2.txt','dfd_rw_140_129_test_input_v2.txt'};
                    {'dfd_rw_140_130_train_input_v2.txt','dfd_rw_140_130_test_input_v2.txt'};
                    {'dfd_rw_140_131_train_input_v2.txt','dfd_rw_140_131_test_input_v2.txt'};
                    {'dfd_rw_140_132_train_input_v2.txt','dfd_rw_140_132_test_input_v2.txt'};
                    {'dfd_rw_140_133_train_input_v2.txt','dfd_rw_140_133_test_input_v2.txt'};
                    {'dfd_rw_140_134_train_input_v2.txt','dfd_rw_140_134_test_input_v2.txt'};
                    {'dfd_rw_140_135_train_input_v2.txt','dfd_rw_140_135_test_input_v2.txt'};
                    {'dfd_rw_140_136_train_input_v2.txt','dfd_rw_140_136_test_input_v2.txt'};
                    {'dfd_rw_140_137_train_input_v2.txt','dfd_rw_140_137_test_input_v2.txt'};
                    {'dfd_rw_140_138_train_input_v2.txt','dfd_rw_140_138_test_input_v2.txt'};
                    {'dfd_rw_140_139_train_input_v2.txt','dfd_rw_140_139_test_input_v2.txt'};
                    {'dfd_rw_139_129_train_input_v2.txt','dfd_rw_139_129_test_input_v2.txt'};
                    {'dfd_rw_139_130_train_input_v2.txt','dfd_rw_139_130_test_input_v2.txt'};
                    {'dfd_rw_139_131_train_input_v2.txt','dfd_rw_139_131_test_input_v2.txt'};
                    {'dfd_rw_139_132_train_input_v2.txt','dfd_rw_139_132_test_input_v2.txt'};
                    {'dfd_rw_139_133_train_input_v2.txt','dfd_rw_139_133_test_input_v2.txt'};
                    {'dfd_rw_139_134_train_input_v2.txt','dfd_rw_139_134_test_input_v2.txt'};
                    {'dfd_rw_139_135_train_input_v2.txt','dfd_rw_139_135_test_input_v2.txt'};
                    {'dfd_rw_139_136_train_input_v2.txt','dfd_rw_139_136_test_input_v2.txt'};
                    {'dfd_rw_139_137_train_input_v2.txt','dfd_rw_139_137_test_input_v2.txt'};
                    {'dfd_rw_139_138_train_input_v2.txt','dfd_rw_139_138_test_input_v2.txt'};
                    {'dfd_rw_138_129_train_input_v2.txt','dfd_rw_138_129_test_input_v2.txt'};
                    {'dfd_rw_138_130_train_input_v2.txt','dfd_rw_138_130_test_input_v2.txt'};
                    {'dfd_rw_138_131_train_input_v2.txt','dfd_rw_138_131_test_input_v2.txt'};
                    {'dfd_rw_138_132_train_input_v2.txt','dfd_rw_138_132_test_input_v2.txt'};
                    {'dfd_rw_138_133_train_input_v2.txt','dfd_rw_138_133_test_input_v2.txt'};
                    {'dfd_rw_138_134_train_input_v2.txt','dfd_rw_138_134_test_input_v2.txt'};
                    {'dfd_rw_138_135_train_input_v2.txt','dfd_rw_138_135_test_input_v2.txt'};
                    {'dfd_rw_138_136_train_input_v2.txt','dfd_rw_138_136_test_input_v2.txt'};
                    {'dfd_rw_138_137_train_input_v2.txt','dfd_rw_138_137_test_input_v2.txt'};
                    {'dfd_rw_137_129_train_input_v2.txt','dfd_rw_137_129_test_input_v2.txt'};
                    {'dfd_rw_137_130_train_input_v2.txt','dfd_rw_137_130_test_input_v2.txt'};
                    {'dfd_rw_137_131_train_input_v2.txt','dfd_rw_137_131_test_input_v2.txt'};
                    {'dfd_rw_137_132_train_input_v2.txt','dfd_rw_137_132_test_input_v2.txt'};
                    {'dfd_rw_137_133_train_input_v2.txt','dfd_rw_137_133_test_input_v2.txt'};
                    {'dfd_rw_137_134_train_input_v2.txt','dfd_rw_137_134_test_input_v2.txt'};
                    {'dfd_rw_137_135_train_input_v2.txt','dfd_rw_137_135_test_input_v2.txt'};
                    {'dfd_rw_137_136_train_input_v2.txt','dfd_rw_137_136_test_input_v2.txt'};
                    {'dfd_rw_136_129_train_input_v2.txt','dfd_rw_136_129_test_input_v2.txt'};
                    {'dfd_rw_136_130_train_input_v2.txt','dfd_rw_136_130_test_input_v2.txt'};
                    {'dfd_rw_136_131_train_input_v2.txt','dfd_rw_136_131_test_input_v2.txt'};
                    {'dfd_rw_136_132_train_input_v2.txt','dfd_rw_136_132_test_input_v2.txt'};
                    {'dfd_rw_136_133_train_input_v2.txt','dfd_rw_136_133_test_input_v2.txt'};
                    {'dfd_rw_136_134_train_input_v2.txt','dfd_rw_136_134_test_input_v2.txt'};
                    {'dfd_rw_136_135_train_input_v2.txt','dfd_rw_136_135_test_input_v2.txt'};
                    {'dfd_rw_135_129_train_input_v2.txt','dfd_rw_135_129_test_input_v2.txt'};
                    {'dfd_rw_135_130_train_input_v2.txt','dfd_rw_135_130_test_input_v2.txt'};
                    {'dfd_rw_135_131_train_input_v2.txt','dfd_rw_135_131_test_input_v2.txt'};
                    {'dfd_rw_135_132_train_input_v2.txt','dfd_rw_135_132_test_input_v2.txt'};
                    {'dfd_rw_135_133_train_input_v2.txt','dfd_rw_135_133_test_input_v2.txt'};
                    {'dfd_rw_135_134_train_input_v2.txt','dfd_rw_135_134_test_input_v2.txt'};
                    {'dfd_rw_134_129_train_input_v2.txt','dfd_rw_134_129_test_input_v2.txt'};
                    {'dfd_rw_134_130_train_input_v2.txt','dfd_rw_134_130_test_input_v2.txt'};
                    {'dfd_rw_134_131_train_input_v2.txt','dfd_rw_134_131_test_input_v2.txt'};
                    {'dfd_rw_134_132_train_input_v2.txt','dfd_rw_134_132_test_input_v2.txt'};
                    {'dfd_rw_134_133_train_input_v2.txt','dfd_rw_134_133_test_input_v2.txt'};
                    {'dfd_rw_133_129_train_input_v2.txt','dfd_rw_133_129_test_input_v2.txt'};
                    {'dfd_rw_133_130_train_input_v2.txt','dfd_rw_133_130_test_input_v2.txt'};
                    {'dfd_rw_133_131_train_input_v2.txt','dfd_rw_133_131_test_input_v2.txt'};
                    {'dfd_rw_133_132_train_input_v2.txt','dfd_rw_133_132_test_input_v2.txt'};
                    {'dfd_rw_132_129_train_input_v2.txt','dfd_rw_132_129_test_input_v2.txt'};
                    {'dfd_rw_132_130_train_input_v2.txt','dfd_rw_132_130_test_input_v2.txt'};
                    {'dfd_rw_132_131_train_input_v2.txt','dfd_rw_132_131_test_input_v2.txt'};
                    {'dfd_rw_131_129_train_input_v2.txt','dfd_rw_131_129_test_input_v2.txt'};
                    {'dfd_rw_131_130_train_input_v2.txt','dfd_rw_131_130_test_input_v2.txt'};
                    {'dfd_rw_130_129_train_input_v2.txt','dfd_rw_130_129_test_input_v2.txt'};
                    {'dfd_rw_dn1_train_input_v2.txt','dfd_rw_dn1_test_input_v2.txt'};
                    {'dfd_rw_dn2_train_input_v2.txt','dfd_rw_dn2_test_input_v2.txt'};
                    {'dfd_rw_dn3_train_input_v2.txt','dfd_rw_dn3_test_input_v2.txt'};
                    {'dfd_rw_up1_train_input_v2.txt','dfd_rw_up1_test_input_v2.txt'};
                    {'dfd_rw_up2_train_input_v2.txt','dfd_rw_up2_test_input_v2.txt'};
                    {'dfd_rw_up3_train_input_v2.txt','dfd_rw_up3_test_input_v2.txt'};
                    {'dfd_rw_up4_train_input_v2.txt','dfd_rw_up4_test_input_v2.txt'}};
    
    file_id = fopen(filename,'w');
    
    fprintf(file_id, '// Auto generated c++ header file for PSO testing\n');
    fprintf(file_id, '// Iteration Number: %d\n', pso_member.iteration);
    fprintf(file_id, '// Population Number: %d\n\n', pso_member.number);

    % print the net headers
    fprintf(file_id, '#ifndef NET_DEFINITION_H\n');
    fprintf(file_id, '#define NET_DEFINITION_H\n\n');
    fprintf(file_id, '#include <cstdint>\n');
    fprintf(file_id, '#include <string>\n\n');
    
    fprintf(file_id, '#include "dlib/dnn.h"\n');
    fprintf(file_id, '#include "dlib/dnn/core.h"\n\n');
    
    fprintf(file_id, '#include "dlib_elu.h"\n');
    fprintf(file_id, '#include "dlib_srelu.h"\n\n');
    
    fprintf(file_id, 'extern const uint32_t img_depth = %d;\n',input_type(pso_member.input,1));
    fprintf(file_id, 'extern const uint32_t secondary = %d;\n',input_type(pso_member.input,2));
    fprintf(file_id, 'extern const std::string train_inputfile = "%s";\n',input_files{pso_member.input_file}{1});
    fprintf(file_id, 'extern const std::string test_inputfile = "%s";\n',input_files{pso_member.input_file}{2});
    fprintf(file_id, 'extern const std::string version = "%s";\n', strcat('pso_',num2str(pso_member.number,'%02d_'),num2str(pso_member.iteration,'%02d_')));
    %fprintf(file_id, 'extern const std::vector<std::pair<uint64_t, uint64_t>> crop_sizes = {{%d, %d}, {30, 140}};\n',pso_member.crop_size*18,pso_member.crop_size*6);
    fprintf(file_id, 'extern const std::vector<std::pair<uint64_t, uint64_t>> crop_sizes = {{%d, %d}, {368, 400}};\n',pso_member.crop_size*crop_scale(1),pso_member.crop_size*crop_scale(2));
    fprintf(file_id, 'extern const uint64_t num_crops = %d;\n',crop_num(pso_member.crop_size));
    fprintf(file_id, 'extern const std::pair<uint32_t, uint32_t> crop_scale(1, 1);\n\n');
    
    fprintf(file_id, '//-----------------------------------------------------------------\n\n');
    fprintf(file_id, 'typedef struct{\n    uint32_t iteration;\n    uint32_t pop_num;\n} pso_struct;\n');
    fprintf(file_id, 'pso_struct pso_info = {%d, %d};\n\n', pso_member.iteration, pso_member.number);

    fprintf(file_id, '//-----------------------------------------------------------------\n\n');
    fprintf(file_id, 'template <long num_filters, typename SUBNET> using con2d = dlib::con<num_filters, 2, 2, 2, 2, SUBNET>;\n');
    fprintf(file_id, 'template <long num_filters, typename SUBNET> using con33d33 = dlib::con<num_filters, 3, 3, 3, 3, SUBNET>;\n');
    fprintf(file_id, 'template <long num_filters, typename SUBNET> using con32d32 = dlib::con<num_filters, 3, 2, 3, 2, SUBNET>;\n');
    fprintf(file_id, 'template <long num_filters, typename SUBNET> using con21d21 = dlib::con<num_filters, 2, 1, 2, 1, SUBNET>;\n\n');

    fprintf(file_id, 'template <long num_filters, typename SUBNET> using cont2u = dlib::cont<num_filters, 2, 2, 2, 2, SUBNET>;\n\n');

    fprintf(file_id, 'template <typename SUBNET> using DTO_0 = dlib::add_tag_layer<200, SUBNET>;\n');
    fprintf(file_id, 'template <typename SUBNET> using DTI_0 = dlib::add_tag_layer<201, SUBNET>;\n');
    fprintf(file_id, 'template <typename SUBNET> using DTO_1 = dlib::add_tag_layer<202, SUBNET>;\n');
    fprintf(file_id, 'template <typename SUBNET> using DTI_1 = dlib::add_tag_layer<203, SUBNET>;\n');
    fprintf(file_id, 'template <typename SUBNET> using DTO_2 = dlib::add_tag_layer<204, SUBNET>;\n');
    fprintf(file_id, 'template <typename SUBNET> using DTI_2 = dlib::add_tag_layer<205, SUBNET>;\n\n');

    fprintf(file_id, '//-----------------------------------------------------------------\n\n');
    fprintf(file_id, 'using dfd_net_type = dlib::loss_multiclass_log_per_pixel<\n');

    for idx=1:numel(net_description.net_structure)

        if(strcmp(net_description.net_structure{idx}, 'con'))
            fprintf(file_id, '    dlib::con<%d, %d, %d, 1, 1, \n', pso_member.con(pso_member.con_map(con_idx),1),2*pso_member.con(pso_member.con_map(con_idx),2)+1,2*pso_member.con(pso_member.con_map(con_idx),3)+1); 
            con_idx = con_idx + 1;
        elseif(strcmp(net_description.net_structure{idx}, 'act'))
            fprintf(file_id, '    %s', net_description.activations{pso_member.act(act_idx)});
            act_idx = act_idx + 1;
        elseif(strcmp(net_description.net_structure{idx}, 'bn'))
            if(pso_member.bn(bn_idx) == 1)
                fprintf(file_id, '    dlib::bn_con<');
            else
                end_cnt = end_cnt + 1;
            end
            bn_idx = bn_idx + 1;
        elseif(strcmp(net_description.net_structure{idx}, 'cond'))
            %fprintf(file_id, '    %s%d,', net_description.cond{cond_order(pso_member.cond,cond_idx)}, pso_member.con(con_idx,1)); 
            fprintf(file_id, '    %s%d,', net_description.cond{cond_idx}, pso_member.con(pso_member.con_map(con_idx),1)); 
            cond_idx = cond_idx + 1;
            con_idx = con_idx + 1;
        elseif(strcmp(net_description.net_structure{idx}, 'cont'))
            fprintf(file_id, '    dlib::cont<%d, %d, %d, 2, 2, \n', pso_member.con(pso_member.con_map(con_idx),1),2*pso_member.con(pso_member.con_map(con_idx),2)+1,2*pso_member.con(pso_member.con_map(con_idx),3)+1);
            con_idx = con_idx + 1;
        elseif(strcmp(net_description.net_structure{idx}, 'cont2u'))
            fprintf(file_id, '    dlib::cont<%d, 2, 2, 2, 2, \n', pso_member.con(pso_member.con_map(con_idx),1));
            con_idx = con_idx + 1;
        elseif(strcmp(net_description.net_structure{idx}, 'concat'))
            fprintf(file_id, '    %s\n', net_description.concats{concat_idx}); 
            concat_idx = concat_idx + 1;        
        elseif(strcmp(net_description.net_structure{idx}, 'add_prev1'))
            fprintf(file_id, '    dlib::add_prev1<');
        elseif(strcmp(net_description.net_structure{idx}, 'tag1'))
            fprintf(file_id, '    dlib::tag1<');
        elseif(strcmp(net_description.net_structure{idx}, 'tags'))
            fprintf(file_id, '    %s',net_description.tags{tag_idx});
            tag_idx = tag_idx + 1;
        elseif(strcmp(net_description.net_structure{idx}, 'input'))
            fprintf(file_id, '    dlib::input<std::array<dlib::matrix<uint16_t>, img_depth>>\n');         
        end    

    end
    
    e = repmat('>',1,idx-end_cnt+1);
    fprintf(file_id, '    %s;\n\n',e);

    fprintf(file_id, '//-----------------------------------------------------------------\n\n');

    fprintf(file_id, 'inline std::ostream& operator<< (std::ostream& out, const pso_struct& item)\n{\n');
    fprintf(file_id, '    out << "------------------------------------------------------------------" << std::endl;\n');
    fprintf(file_id, '    out << "PSO Info: " << std::endl;\n');
    fprintf(file_id, '    out << "  Iteration: " << item.iteration << std::endl;\n');
    fprintf(file_id, '    out << "  Population Number: " << item.pop_num << std::endl;\n');
    fprintf(file_id, '    out << "------------------------------------------------------------------" << std::endl;\n');
    fprintf(file_id, '    return out;\n}\n\n');

    fprintf(file_id, '#endif \n');
    
    fclose(file_id); 

end