#include <memory>
#include <restbed>
#include "nexson.h"
#include "otc/error.h"
#include "ws/tolwsadaptors.h"

using std::string;
using std::pair;
using nlohmann::json;
using std::make_unique;
using std::make_shared;
using namespace restbed;

namespace otc {

// https://github.com/OpenTreeOfLife/reference-taxonomy/blob/master/org/opentreeoflife/server/Services.java#L285

// This should probably be a parameter
static string studybase = "https://api.opentreeoflife.org/v3/study/";
// reference-taxonomy also caches a single study.  We could add a use_cache argument to do the same.
json get_phylesystem_study(const string& study_id)
{
    using namespace restbed;

    LOG(WARNING)<<"getting phylesystem study '"<<study_id<<"'";

    // 0. Construct session
    auto session = make_shared< Session > ("");

    // 1. Construct request
    Uri uri( studybase + study_id + "?output_nexml2json=1.2.1");
    auto request = make_shared< Request >( uri );
    request->set_header( "Accept", "*/*" );
    request->set_header( "Host", "api.opentreeoflife.org" );
    request->set_header( "User-Agent", "otc-tol-ws" );
//    request->set_query_parameter( "output_nexml2json", "1.2.1" );
//    request->set_query_parameter( "auth_token", "ANONYMOUS" );

    // 2. Send request
    LOG(WARNING)<<"reading uri "<<uri.to_string();
//    session->send(request, 
    auto response = Http::sync( request );


    LOG(WARNING)<<"Status Code:    "<<response->get_status_code()<<"\n";
    LOG(WARNING)<<"Status Message: "<<response->get_status_message().data()<<"\n";
    LOG(WARNING)<<"HTTP Version:   "<<response->get_version()<<"\n";
    LOG(WARNING)<<"HTTP Protocol:  "<<response->get_protocol().data()<<"\n";

    for ( const auto header : response->get_headers( ) )
        LOG(WARNING)<<"Header '"<<header.first.data()<<"' > '"<<header.second.data()<<"'";

    if (response->get_status_code() != 200)
	throw OTCBadRequest()<<"GET '"<<uri.to_string()<<"' yielded "<<response->get_status_code();

    // 3. Read request
    if (response->has_header("Transfer-Encoding"))
    {
        // https://github.com/Corvusoft/restbed/blob/master/example/transfer_encoding_request/source/example.cpp
	LOG(WARNING)<<"got HERE 2a";
	Http::fetch("\r\n", response);
	const auto & data = response->get_body();
	if (not data.empty())
	{
	    const string length (data.begin(), data.end());
	    if (length != "0\r\n")
	    {
		const auto chunk_size = stoul(length, nullptr, 16) + strlen("\r\n");
		LOG(WARNING)<<"chunk size = "<<chunk_size;
	    }
	}
	LOG(WARNING)<<"got HERE 2a.  response has length "<<response->get_body().size();
    }
    else
    {
	auto length = response->get_header( "Content-Length", 0 );
	LOG(WARNING)<<"got HERE 2b.  Content-Length = "<<length;

	Http::fetch(length, response);
	LOG(WARNING)<<"got HERE 2b.  response has length "<<response->get_body().size();
    }

    // 4. Convert quest to JSON
    auto j = parse_body( response->get_body() );
    if (not j)
	throw OTCBadRequest()<<"Could not parse JSON for study "<<study_id;

    if (not j->count("data"))
	throw OTCBadRequest()<<"No 'data' property in response to GET "<<uri.to_string();

    j = (*j)["data"];

    if (not j->count("nexml"))
	throw OTCBadRequest()<<"No 'nexml' property in JSON data blob from "<<uri.to_string();

    return *j;
}

// See extract_tree_nexson in peyotl/nexson_syntax/__init__.py
pair<json,json> extract_tree_nexson(const json& nexson, const string& treeid)
{
    if (not nexson.count("nexml"))
	throw OTCError()<<"No 'nexml' element in json blob";
    auto nexml_el = nexson["nexml"];

    if (not nexml_el.count("treesById"))
	throw OTCError()<<"No 'treesById' element in nexml element";
    json tree_groups = nexml_el["treesById"];

    for(auto& tree_group: tree_groups)
    {
	auto trees = tree_group["treeById"];
	if (not trees.count(treeid)) continue;

	auto tree = trees[treeid];
	auto otu_groups = nexml_el["otusById"];
	string ogi = tree_group["@otus"];
	auto otu_group = otu_groups[ogi]["otuById"];
	return {tree, otu_group};
    }
    throw OTCError()<<"No tree '"<<treeid<<"' found in study";
}

} //namespace otc
