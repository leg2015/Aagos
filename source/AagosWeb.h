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

  UI::Button run_toggle_but;
  UI::Button run_step_but;

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

  /// Advance application by one frame.
  void DoFrame();
};

void AagosWebInterface::SetupInterface() {
  world_div << "Hello world";

  // ---- Setup interface control buttons ----
  // Run button setup.
  // Manually create run/pause button to have full control over button's callback.
  // NOTE - here's an example where Empirical's web tools are overly restrictive.
  run_toggle_but = UI::Button([this]() {
    ToggleActive();
    run_toggle_but.SetLabel(active ? "Stop" : "Run"); // TODO - use icons
    active ? run_step_but.SetDisabled(true) : run_step_but.SetDisabled(false);
  }, "Run", "run-toggle-button");
  run_toggle_but.SetAttr("class", "btn btn-primary m-1 ");

  // Step button setup.
  run_step_but = GetStepButton("run-step-button");
  run_step_but.SetAttr("class", "btn btn-primary m-1");

  // Add buttons to
  control_div << run_step_but;
  control_div << run_toggle_but;
}

void AagosWebInterface::DoFrame() {
  std::cout << "Frame!" << std::endl;
}

#endif