#include <cerrno>
#include <cstdarg>
#include <cstddef>
#include <cstdlib>
#include <string>

extern "C" {
#include <libxml/xmlerror.h>
}

#include <arbornml/nmlexcept.hpp>

#include "xmlwrap.hpp"

namespace arbnml {

bool parse_attr(std::string& out, const char* content) {
    out = content;
    return true;
}

bool parse_attr(long long& out, const char* content) {
    char* end = 0;
    out = std::strtoll(content, &end, 10);
    if (end && !*end) return true;
    return false;
}

bool parse_attr(unsigned long long& out, const char* content) {
    if (std::strtol(content, nullptr, 10)>=0) {
        errno = 0;
        char* end = 0;
        out = std::strtoull(content, &end, 10);
        if (!errno && end && !*end) return true;
    }
    return false;
}

bool parse_attr(double& out, const char* content) {
    char* end = 0;
    out = std::strtod(content, &end);
    return end && !*end;
}

void throw_on_xml_generic_error(void *, const char* msg, ...) {
    va_list va, vb;
    va_start(va, msg);
    va_copy(vb, va);

    int r = vsnprintf(nullptr, 0, msg, va);
    va_end(va);

    std::string err(r+1, '\0');
    vsnprintf(&err[0], err.size(), msg, vb);
    va_end(vb);

    throw ::arbnml::xml_error(err);
}

void throw_on_xml_structured_error(void *ctx, xmlErrorPtr errp) {
    if (errp->level!=1) { // ignore warnings!
        std::string msg(errp->message);
        if (!msg.empty() && msg.back()=='\n') msg.pop_back();
        throw ::arbnml::xml_error(msg, errp->line);
    }
}

xml_error_scope::xml_error_scope() {
    generic_handler_ = xmlGenericError;
    generic_context_ = xmlGenericErrorContext;

    structured_handler_ = xmlStructuredError;
    structured_context_ = xmlStructuredErrorContext;

    xmlSetGenericErrorFunc(nullptr, &throw_on_xml_generic_error);
    xmlSetStructuredErrorFunc((void*)this, &throw_on_xml_structured_error);
}

xml_error_scope::~xml_error_scope() {
    xmlGenericError = generic_handler_;
    xmlGenericErrorContext = generic_context_;

    xmlStructuredError = structured_handler_;
    xmlStructuredErrorContext = structured_context_;
}

} // namespace arbnml
