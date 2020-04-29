#ifndef AAGOS_WEB_H
#define AAGOS_WEB_H

#include "web/web.h"

#include "AagosWorld.h"
#include "AagosConfig.h"

namespace UI = emp::web;

double GetHTMLElementWidthByID(const std::string & id) {
  return EM_ASM_DOUBLE({
      var id = UTF8ToString($0);
      return $('#' + id).width();
    }, id.c_str());
}

double GetHTMLElementHeightByID(const std::string & id) {
  return EM_ASM_DOUBLE({
      var id = UTF8ToString($0);
      return $('#' + id).height();
    }, id.c_str());
}

class AagosWebInterface : public UI::Animate, public AagosWorld {
public:
  using config_t = AagosConfig;

protected:
public:
  AagosWebInterface(config_t & cfg) : AagosWorld(cfg) {
    std::cout << "AagowWebInterface constructor!" << std::endl;
  }
};

#endif