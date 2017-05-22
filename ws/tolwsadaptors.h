#ifndef OTC_TOLWS_ADAPTORS_H
#define OTC_TOLWS_ADAPTORS_H
#include <restbed>
#include "ws/tolws.h"
#include "ws/parallelreadserialwrite.h"
#include "otc/otc_base_includes.h"

template<typename T>
bool extract_from_request(const nlohmann::json & j,
                          std::string opt_name,
                          T & setting,
                          std::string & response,
                          int & status_code);

inline bool set_type_expectations_error(const std::string & opt_name,
                                        const char * type_name_with_article,
                                        std::string & response,
                                        int & status_code) {
    response = "Expecting ";
    response += opt_name;
    response += " to be ";
    response += type_name_with_article;
    response += ".\n";
    status_code = 400;
    return false;
}

template<>
inline bool extract_from_request(const nlohmann::json & j,
                          std::string opt_name,
                          bool & setting,
                          std::string & response,
                          int & status_code) {
    auto opt = j.find(opt_name);
    if (opt != j.end()) {
        if (opt->is_boolean()) {
            setting = opt->get<bool>();
            return true;
        }
        return set_type_expectations_error(opt_name, "a boolean", response, status_code);
    }
    return false;
}

template<>
inline bool extract_from_request(const nlohmann::json & j,
                          std::string opt_name,
                          std::string & setting,
                          std::string & response,
                          int & status_code) {
    auto opt = j.find(opt_name);
    if (opt != j.end()) {
        if (opt->is_string()) {
            setting = opt->get<std::string>();
            return true;
        }
        return set_type_expectations_error(opt_name, "a string", response, status_code);
    }
    return false;
}

template<>
inline bool extract_from_request(const nlohmann::json & j,
                          std::string opt_name,
                          int & setting,
                          std::string & response,
                          int & status_code) {
    auto opt = j.find(opt_name);
    if (opt != j.end()) {
        if (opt->is_number()) {
            long ls = opt->get<long>();
            if (ls >= std::numeric_limits<int>::max()) {
              return set_type_expectations_error(opt_name, "a non-huge integer", response, status_code);
            }
            setting = static_cast<int>(ls);
            return true;
        }
        return set_type_expectations_error(opt_name, "an integer", response, status_code);
    }
    return false;
}

template<>
inline bool extract_from_request(const nlohmann::json & j,
                          std::string opt_name,
                          long & setting,
                          std::string & response,
                          int & status_code) {
    auto opt = j.find(opt_name);
    if (opt != j.end()) {
        if (opt->is_number()) {
            setting = opt->get<long>();
            return true;
        }
        return set_type_expectations_error(opt_name, "an integer", response, status_code);
    }
    return false;
}

template<>
inline bool extract_from_request(const nlohmann::json & j,
                          std::string opt_name,
                          std::vector<std::string> & setting,
                          std::string & response,
                          int & status_code) {
    auto opt = j.find(opt_name);
    if (opt != j.end()) {
        if (opt->is_array()) {
            for (nlohmann::json::const_iterator sIt = opt->begin(); sIt != opt->end(); ++sIt) {
                if (sIt->is_string()) {
                    setting.push_back(sIt->get<std::string>());
                } else {
                    return set_type_expectations_error(opt_name, "an array of strings", response, status_code);
                }   
            }
            return true;
        }
        return set_type_expectations_error(opt_name, "an array of strings", response, status_code);
    }
    return false;
}

template<>
inline bool extract_from_request(const nlohmann::json & j,
                                 std::string opt_name,
                                 otc::OttIdSet & setting,
                                 std::string & response,
                                 int & status_code) {
    auto opt = j.find(opt_name);
    if (opt != j.end()) {
        if (opt->is_array()) {
            for (nlohmann::json::const_iterator sIt = opt->begin(); sIt != opt->end(); ++sIt) {
                if (sIt->is_number()) {
                    setting.insert(sIt->get<otc::OttId>());
                } else {
                    return set_type_expectations_error(opt_name, "an array of integers", response, status_code);
                }   
            }
            return true;
        }
        return set_type_expectations_error(opt_name, "an array of integers", response, status_code);
    }
    return false;
}



inline otc::vec_src_node_ids extract_node_id_vec(otc::TreesToServe & tts,
                                                 const nlohmann::json & sbv
#                                                if defined(JOINT_MAPPING_VEC)
                                                   , otc::SourceEdgeMappingType semt
#                                                endif
                                                 ) {
#   if defined(JOINT_MAPPING_VEC)
       using lel_t = otc::semt_ind_t;
#   else
       using lel_t = std::uint32_t;
#   endif
    std::list<lel_t> lsni;
    for (nlohmann::json::const_iterator jit = sbv.begin(); jit != sbv.end(); ++jit) {
        const std::string * kp = tts.get_stored_string(jit.key());
        const auto & v = jit.value();
        for (nlohmann::json::const_iterator vit = v.begin(); vit != v.end(); ++vit) {
            const std::string * vp = tts.get_stored_string(*vit);
            const auto sni_ind = tts.get_source_node_id_index(otc::src_node_id(kp, vp));
#           if defined(JOINT_MAPPING_VEC)
                lsni.push_back(lel_t(semt, sni_ind));
#           else
                lsni.push_back(sni_ind);
#           endif
        } 
    }
    return otc::vec_src_node_ids(lsni.begin(), lsni.end());
}


inline const nlohmann::json & extract_obj(const nlohmann::json &j, const char * field) {
    auto dc_el = j.find(field);
    if (dc_el == j.end()) {
        throw otc::OTCError() << "Missing \"" << field << "\" field.\n";
    }
    if (dc_el->is_object()) {
        nlohmann::json & k = const_cast<nlohmann::json &>(j);
        return k[field];
    }
    throw otc::OTCError() << "Expected \"" << field << "\" field to be a string.\n";
}


///////////////////////
// handlers that are registered as callback
void mrca_method_handler(const std::shared_ptr<restbed::Session> session);
void induced_subtree_method_handler(const std::shared_ptr<restbed::Session> session);
void tax_about_method_handler(const std::shared_ptr<restbed::Session> session);
void taxon_info_method_handler(const std::shared_ptr<restbed::Session> session);
int run_server(const boost::program_options::variables_map & args);
boost::program_options::variables_map parse_cmd_line(int argc, char* argv[]);

#endif
