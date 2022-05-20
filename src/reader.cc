#include "reader.hh"

#include <iostream>

#include "utils.hh"
#include "toml.hpp"

// Suppress optimization to speed up compilation
#if defined(__GNUC__) && __GNUC__ >=4
#pragma GCC optimize "-O"
#endif

// Suppress spurious warning triggered by toml11 for gcc 4.8.5
#if defined(__GNUC__) && __GNUC__ < 5
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#include <fstream>

// Re-enable warnings
#if defined(__GNUC__) && __GNUC__ < 5
#pragma GCC diagnostic warning "-Wunused-parameter"
#endif

// Wrapper class for TOML type
class Node {
public:
	typedef toml::basic_value<toml::discard_comments,
														std::unordered_map, std::vector> value_type;
	explicit Node(const value_type& v) : v(v) {}
	value_type v;
};

using namespace std;

InputNode::InputNode(
	const Node n,
	const std::string& label,
	const InputNode* parent
	) :
	n_(std::shared_ptr<Node>(new Node{n})), label_(label), parent_(parent)
{
}

size_t InputNode::size() const
{
	return n().v.size();
}

const Node& InputNode::n() const
{
	return *n_;
}

bool InputNode::contains(const std::string& name) const
{
	return n().v.contains(name);
}

// Use index to augment label
InputNode InputNode::operator[](
	unsigned int index) const
{
	std::ostringstream oss;
	oss << '[' << index << ']';
	return InputNode(Node(toml::find(n().v,index)),oss.str(),this);
}

InputNode InputNode::operator[](
	const std::string& s) const
{
	return InputNode(Node(toml::find(n().v,s)),s,this);
}

void InputNode::set_used()
{
	for(const auto &p : n().v.as_table()){
		UsedTomlKey use; use.name = p.first; use.used = false;
		used.push_back(use);
	}
}

void InputNode::add_used(const string &name)
{
	auto i = 0u; while(i < used.size() && used[i].name != name) i++;
	if(i < used.size()) used[i].used = true;
}

void InputNode::check_used(string name)
{
	vector <string> not_used;	
	for(const auto &use : used){
		if(use.used == false) not_used.push_back(use.name);
	}

	if(not_used.size() > 0){
		int num;
		MPI_Comm_rank(MPI_COMM_WORLD,&num);
		if(num == 0){
			cout << "WARNING! In '" << name << "' ";
			switch(not_used.size()){
				case 1: cout << "the quantity '" << not_used[0] << "' is not used."; break;
				default: 
					cout << "the quantities ";
					for(auto i = 0u; i < not_used.size(); i++){
						if(i != 0){
							if(i < not_used.size()-1) cout << ", "; else cout << " and ";
						}
						cout << "'" << not_used[i] << "'";		
					}
					cout << " are not used.";
					break;
			}
			cout << endl << endl;
		}
	}
}

std::string InputNode::stringfield_unchecked(const std::string& name)
{
	add_used(name);
	auto st = toml::find<std::string>(n().v,name);
	strip(st);
	auto i =0u; while(i < st.length() && st.substr(i,1) != ",") i++;
	if(i < st.length()) emsgroot("The expression '"+st+"' in '"+name+"' must not contain a comma");
	return st;
}

std::vector <double> InputNode::find_vector(const std::string& name)
{
	add_used(name);
  return toml::find<std::vector<double>>(n().v,name);	
}

std::string InputNode::stringfield(string name, string em)
{
	return stringfield(name.c_str(),em);
}

std::string InputNode::stringfield(const char *name, string em)
{
	if(!contains(name)){
		if(em == "") return "";
		else emsgroot(em+" a value for '"+name+"' must be specified");
	}
	return stringfield_unchecked(name);
}

std::string InputNode::stringfield_allowcomma(string name, string em)
{
	return stringfield_allowcomma(name.c_str(),em);
}

std::string InputNode::stringfield_allowcomma(const char *name, string em)
{
	if(!contains(name)){
		if(em == "") return "";
		else emsgroot(em+" a value for '"+name+"' must be specified");
	}
	add_used(name);
	auto st = toml::find<std::string>(n().v,name);
	strip(st);
	return st;
}

double InputNode::numberfield_unchecked(const std::string& name)
{
	const auto value = toml::find(n().v,name);
	if(value.is_floating()) return value.as_floating();
	else return value.as_integer();
}

double InputNode::numberfield(const char *name, string em)
{
	if(!contains(name)){
		if(em == "") return UNSET;
		else emsgroot(em+" a value for '"+name+"' must be specified");
	}
	return numberfield_unchecked(name);
}

int InputNode::intfield_unchecked(const std::string& name)
{
	add_used(name);
	return toml::find<int>(n().v,name);
}

/// Gets a list of all the keys
vector<string> InputNode::keys( ) const
{
	vector<string> keys;
	for(const auto &p : n().v.as_table())
	{
		keys.push_back(p.first);
	}
	return keys;
}



/// Factory to generate InputNode from names file
InputNode parsefile(const std::string& inputfilename)
{
	Node n{toml::parse(inputfilename)};
	// Allow using values from another TOML file as a base for this one. TODO:
	// make this into functions so you can do this recursively.
	if (n.v.contains("baseinputfile")) {
		const string basefile = toml::find<std::string>(n.v,"baseinputfile");
		// TODO: make the filename relative to the original TOML file
		decltype(toml::parse(basefile)) basetomlddata = toml::parse(basefile);
		
		for(const auto &p : basetomlddata.as_table())
		{
			if (!n.v.contains(p.first)) {
				n.v[p.first] = p.second;
			}
		}
	}
	return InputNode(n,"",NULL);
}
