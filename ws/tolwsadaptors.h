#ifndef OTC_TOLWS_ADAPTORS_H
#define OTC_TOLWS_ADAPTORS_H
#include <restbed>
#include "ws/tolws.h"
template<typename T>
bool extract_from_request(const nlohmann::json & j,
                          std::string opt_name,
                          T & setting,
                          std::string & response,
                          int & status_code);

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
        response = "Expecting ";
        response += opt_name;
        response += " to be a boolean.\n";
        status_code = 400;
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
        response = "Expecting ";
        response += opt_name;
        response += " to be a string.\n";
        status_code = 400;
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
            setting = opt->get<int>();
            return true;
        }
        response = "Expecting ";
        response += opt_name;
        response += " to be a number.\n";
        status_code = 400;
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
        response = "Expecting ";
        response += opt_name;
        response += " to be a number.\n";
        status_code = 400;
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
                    response = "Expecting ";
                    response += opt_name;
                    response += " to be an array of strings.\n";
                    status_code = 400;
                    return false;
                }   
            }
            return true;
        }
        response = "Expecting ";
        response += opt_name;
        response += " to be an array of strings.\n";
        status_code = 400;
    }
    return false;
}


inline otc::vec_src_node_ids extract_node_id_vec(otc::TreesToServe & tts,
                                                 const nlohmann::json & sbv) {
    std::list<otc::src_node_id> lsni;
    for (nlohmann::json::const_iterator jit = sbv.begin(); jit != sbv.end(); ++jit) {
        const std::string * kp = tts.getStoredString(jit.key());
        const auto & v = jit.value();
        for (nlohmann::json::const_iterator vit = v.begin(); vit != v.end(); ++vit) {
            const std::string * vp = tts.getStoredString(*vit);
            lsni.push_back(otc::src_node_id(kp, vp));
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
void about_method_handler(const std::shared_ptr<restbed::Session> session);
void node_info_method_handler(const std::shared_ptr<restbed::Session> session);
void mrca_method_handler(const std::shared_ptr<restbed::Session> session);
void subtree_method_handler(const std::shared_ptr<restbed::Session> session);
void induced_subtree_method_handler(const std::shared_ptr<restbed::Session> session);
void tax_about_method_handler(const std::shared_ptr<restbed::Session> session);
void taxon_info_method_handler(const std::shared_ptr<restbed::Session> session);
int run_server(const boost::program_options::variables_map & args);
boost::program_options::variables_map parse_cmd_line(int argc, char* argv[]);

#endif