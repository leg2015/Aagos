#ifndef AAGOS_WEB_H
#define AAGOS_WEB_H

/**
 * TODOs
 * - [ ] Allow world to be configured by web user

*/

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

  UI::Document world_div;
  UI::Document control_div;

public:
  AagosWebInterface(config_t & cfg)
    : AagosWorld(cfg),
      world_div("emp_world_view"),
      control_div("emp_controls_view")
  {
    std::cout << "AagowWebInterface constructor!" << std::endl;
    SetupInterface();
    // Setup world.
    Setup();
  }

  /// Configure the interface
  void SetupInterface();
};

void AagosWebInterface::SetupInterface() {
  world_div << "Hello world";
  control_div << "Hello again";

}

#endif