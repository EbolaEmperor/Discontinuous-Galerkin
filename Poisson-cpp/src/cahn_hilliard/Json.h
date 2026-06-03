#ifndef CH_JSON_H
#define CH_JSON_H

// Minimal, dependency-free JSON reader for flat config files.
// Supports the standard JSON grammar plus // line and /* */ block comments
// (so config files can be annotated). Intended for small key->scalar configs;
// nested objects/arrays parse fine but the typed getters below target scalars.

#include <string>
#include <map>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cctype>

namespace cfgjson {

struct Json {
    enum Type { Null, Bool, Number, String, Array, Object };
    Type type = Null;
    bool b = false;
    double num = 0.0;
    std::string str;
    std::vector<Json> arr;
    std::map<std::string, Json> obj;

    bool contains(const std::string& k) const {
        return type == Object && obj.find(k) != obj.end();
    }
    bool hasNumber(const std::string& k) const {
        return contains(k) && obj.at(k).type == Number;
    }
    double getNumber(const std::string& k, double def) const {
        if (!contains(k)) return def;
        const Json& v = obj.at(k);
        if (v.type == Number) return v.num;
        if (v.type == Bool)   return v.b ? 1.0 : 0.0;
        return def;
    }
    int getInt(const std::string& k, int def) const {
        return static_cast<int>(getNumber(k, static_cast<double>(def)));
    }
    long getLong(const std::string& k, long def) const {
        return static_cast<long>(getNumber(k, static_cast<double>(def)));
    }
    bool getBool(const std::string& k, bool def) const {
        if (!contains(k)) return def;
        const Json& v = obj.at(k);
        if (v.type == Bool)   return v.b;
        if (v.type == Number) return v.num != 0.0;
        return def;
    }
    std::string getString(const std::string& k, const std::string& def) const {
        if (!contains(k)) return def;
        const Json& v = obj.at(k);
        return (v.type == String) ? v.str : def;
    }
};

class Parser {
public:
    explicit Parser(const std::string& s) : s_(s), i_(0) {}
    Json parse() {
        Json v = parseValue();
        skipWs();
        if (i_ != s_.size()) fail("trailing characters after JSON value");
        return v;
    }
private:
    const std::string& s_;
    size_t i_;

    [[noreturn]] void fail(const std::string& msg) {
        int line = 1, col = 1;
        for (size_t k = 0; k < i_ && k < s_.size(); ++k) {
            if (s_[k] == '\n') { ++line; col = 1; } else ++col;
        }
        throw std::runtime_error("line " + std::to_string(line) + " col " +
                                 std::to_string(col) + ": " + msg);
    }
    char peek() const { return i_ < s_.size() ? s_[i_] : '\0'; }

    void skipWs() {
        for (;;) {
            while (i_ < s_.size() && std::isspace(static_cast<unsigned char>(s_[i_]))) ++i_;
            if (i_ + 1 < s_.size() && s_[i_] == '/' && s_[i_ + 1] == '/') {
                i_ += 2;
                while (i_ < s_.size() && s_[i_] != '\n') ++i_;
            } else if (i_ + 1 < s_.size() && s_[i_] == '/' && s_[i_ + 1] == '*') {
                i_ += 2;
                while (i_ + 1 < s_.size() && !(s_[i_] == '*' && s_[i_ + 1] == '/')) ++i_;
                if (i_ + 1 < s_.size()) i_ += 2; else fail("unterminated block comment");
            } else {
                break;
            }
        }
    }

    Json parseValue() {
        skipWs();
        char c = peek();
        if (c == '{') return parseObject();
        if (c == '[') return parseArray();
        if (c == '"') { Json v; v.type = Json::String; v.str = parseString(); return v; }
        if (c == '-' || (c >= '0' && c <= '9')) return parseNumber();
        if (c == 't' || c == 'f') return parseLiteralBool();
        if (c == 'n') return parseLiteralNull();
        fail("unexpected character");
    }

    Json parseObject() {
        Json v; v.type = Json::Object;
        ++i_; // consume '{'
        skipWs();
        if (peek() == '}') { ++i_; return v; }
        for (;;) {
            skipWs();
            if (peek() != '"') fail("expected string key");
            std::string key = parseString();
            skipWs();
            if (peek() != ':') fail("expected ':'");
            ++i_;
            v.obj[key] = parseValue();
            skipWs();
            char c = peek();
            if (c == ',') { ++i_; continue; }
            if (c == '}') { ++i_; break; }
            fail("expected ',' or '}'");
        }
        return v;
    }

    Json parseArray() {
        Json v; v.type = Json::Array;
        ++i_; // consume '['
        skipWs();
        if (peek() == ']') { ++i_; return v; }
        for (;;) {
            v.arr.push_back(parseValue());
            skipWs();
            char c = peek();
            if (c == ',') { ++i_; continue; }
            if (c == ']') { ++i_; break; }
            fail("expected ',' or ']'");
        }
        return v;
    }

    std::string parseString() {
        ++i_; // consume opening quote
        std::string out;
        for (;;) {
            if (i_ >= s_.size()) fail("unterminated string");
            char c = s_[i_++];
            if (c == '"') break;
            if (c == '\\') {
                if (i_ >= s_.size()) fail("unterminated escape");
                char e = s_[i_++];
                switch (e) {
                    case '"':  out += '"';  break;
                    case '\\': out += '\\'; break;
                    case '/':  out += '/';  break;
                    case 'b':  out += '\b'; break;
                    case 'f':  out += '\f'; break;
                    case 'n':  out += '\n'; break;
                    case 'r':  out += '\r'; break;
                    case 't':  out += '\t'; break;
                    case 'u': {
                        if (i_ + 4 > s_.size()) fail("bad \\u escape");
                        unsigned code = 0;
                        for (int k = 0; k < 4; ++k) {
                            char h = s_[i_++];
                            code <<= 4;
                            if (h >= '0' && h <= '9') code |= unsigned(h - '0');
                            else if (h >= 'a' && h <= 'f') code |= unsigned(h - 'a' + 10);
                            else if (h >= 'A' && h <= 'F') code |= unsigned(h - 'A' + 10);
                            else fail("bad hex in \\u escape");
                        }
                        if (code < 0x80) {
                            out += static_cast<char>(code);
                        } else if (code < 0x800) {
                            out += static_cast<char>(0xC0 | (code >> 6));
                            out += static_cast<char>(0x80 | (code & 0x3F));
                        } else {
                            out += static_cast<char>(0xE0 | (code >> 12));
                            out += static_cast<char>(0x80 | ((code >> 6) & 0x3F));
                            out += static_cast<char>(0x80 | (code & 0x3F));
                        }
                        break;
                    }
                    default: fail("invalid escape");
                }
            } else {
                out += c;
            }
        }
        return out;
    }

    Json parseNumber() {
        size_t start = i_;
        if (peek() == '-') ++i_;
        while (std::isdigit(static_cast<unsigned char>(peek()))) ++i_;
        if (peek() == '.') { ++i_; while (std::isdigit(static_cast<unsigned char>(peek()))) ++i_; }
        if (peek() == 'e' || peek() == 'E') {
            ++i_;
            if (peek() == '+' || peek() == '-') ++i_;
            while (std::isdigit(static_cast<unsigned char>(peek()))) ++i_;
        }
        Json v; v.type = Json::Number;
        v.num = std::strtod(s_.substr(start, i_ - start).c_str(), nullptr);
        return v;
    }

    Json parseLiteralBool() {
        if (s_.compare(i_, 4, "true") == 0)  { i_ += 4; Json v; v.type = Json::Bool; v.b = true;  return v; }
        if (s_.compare(i_, 5, "false") == 0) { i_ += 5; Json v; v.type = Json::Bool; v.b = false; return v; }
        fail("invalid literal (expected true/false)");
    }
    Json parseLiteralNull() {
        if (s_.compare(i_, 4, "null") == 0) { i_ += 4; Json v; v.type = Json::Null; return v; }
        fail("invalid literal (expected null)");
    }
};

// Parse JSON text; throws std::runtime_error on malformed input.
inline Json parse(const std::string& text) { return Parser(text).parse(); }

// Read a whole file into `text`. Returns false if the file cannot be opened.
inline bool readFile(const std::string& path, std::string& text) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    std::stringstream ss;
    ss << f.rdbuf();
    text = ss.str();
    return true;
}

} // namespace cfgjson

#endif
