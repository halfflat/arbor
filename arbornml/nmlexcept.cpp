#include <string>

#include <arbornml/nmlexcept.hpp>

namespace arbnml {

static std::string fmt_error(const char* prefix, const std::string& err, unsigned line) {
    return prefix + (line==0? err: "line " + std::to_string(line) + ": " + err);
}

xml_error::xml_error(const std::string& xml_error_msg, unsigned line):
    neuroml_exception(fmt_error("xml error: ", xml_error_msg, line)),
    xml_error_msg(xml_error_msg),
    line(line)
{}

no_document::no_document():
    neuroml_exception("no NeuroML document to parse")
{}

parse_error::parse_error(const std::string& error_msg, unsigned line):
    neuroml_exception(fmt_error("parse error: ", error_msg, line)),
    error_msg(error_msg),
    line(line)
{}

bad_segment::bad_segment(unsigned long long segment_id, unsigned line):
    neuroml_exception(
        fmt_error(
            "bad morphology segment: ",
            "segment "+(segment_id+1==0? "unknown": "\""+std::to_string(segment_id)+"\""),
            line)),
    segment_id(segment_id),
    line(line)
{}

} // namespace arbnml
