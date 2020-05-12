#ifndef AAGOS_WEB_D3_H
#define AAGOS_WEB_D3_H

/**
 * TODOs
 * - [ ] Allow world to be configured by web user
*/


#include "web/web.h"
#include "web/Input.h"
#include "web/Element.h"
#include "tools/string_utils.h"

#include <unordered_map>

#include "AagosWorld.h"
#include "AagosConfig.h"
#include "AagosOrg.h"
#include "AagosPopulationVisualization.h"

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
  // static constexpr double BIT_WIDTH
  static constexpr double BIT_HEIGHT = 20;
  static constexpr double MAX_GENE_TARGET_BIT_WIDTH = 10;
  static constexpr double GENE_TARGET_IDENTIFER_HEIGHT = 5;
  static constexpr double INDIV_VERT_MARGIN = 5;

  using config_t = AagosConfig;
  using canvas_draw_fun_t = std::function<void(void)>;
  using input_callback_fun_t = std::function<void(std::string)>;
  using input_checker_fun_t = std::function<bool(std::string)>;

  using org_t = AagosOrg;
  using genome_t = AagosOrg::Genome;
  using phenotype_t = AagosOrg::Phenotype;

protected:

  UI::Document world_div;
  UI::Document control_div;
  UI::Document config_general_div;
  UI::Document config_genetic_arch_div;
  UI::Document config_mutation_div;
  UI::Document config_phase_2_evo_div;

  UI::Button run_toggle_but;
  UI::Button run_step_but;

  std::unordered_map<std::string, UI::Input> config_input_elements;

  size_t draw_frequency=32;

  AagosPopulationVisualization pop_vis;

  void RedrawPopulation(bool update_data=true);
  void RedrawEnvironment();

  void SetupConfigInterface();

  void DisableConfigInputs(bool in_dis=true) {
    for (auto & cfg : config_input_elements) {
      cfg.second.Disabled(in_dis);
    }
  }

  bool CheckInputSize_t(const std::string & in) {
    return (emp::is_digits(in) && in.size() && std::stoi(in) > 0);
  }

public:
  AagosWebInterface(config_t & cfg)
    : AagosWorld(cfg),
      world_div("emp_world_view"),
      control_div("emp_controls_view"),
      config_general_div("emp_view_config_general"),
      config_genetic_arch_div("emp_view_config_genetic_architecture"),
      config_mutation_div("emp_view_config_mutation"),
      config_phase_2_evo_div("emp_view_config_phase_2_evolution"),
      pop_vis("emp-pop-vis")
  {
    std::cout << "AagowWebInterface constructor!" << std::endl;
    SetupInterface();
  }

  /// Configure the interface
  void SetupInterface();

  /// Advance application by one frame.
  void DoFrame();
};

void AagosWebInterface::SetupInterface() {
  // ---- Setup interface control buttons ----
  // Run button setup.
  // Manually create run/pause button to have full control over button's callback.
  // NOTE - here's an example where Empirical's web tools are overly restrictive.
  run_toggle_but = UI::Button([this]() {
    ToggleActive();
    run_toggle_but.SetLabel(active ? "Stop" : "Run"); // TODO - use icons
    run_step_but.SetDisabled(active);
    DisableConfigInputs(active);
  }, "Run", "run-toggle-button");
  run_toggle_but.SetAttr("class", "btn btn-block btn-lg btn-primary");
  std::cout << "setup run_toggle_but" << std::endl;

  // Step button setup.
  run_step_but = GetStepButton("run-step-button");
  run_step_but.SetAttr("class", "btn btn-block btn-lg btn-primary");

  // Add buttons to controls view.
  control_div << UI::Div("button-row").SetAttr("class", "row justify-content-md-center");

  control_div.Div("button-row")
    << UI::Div("step-col").SetAttr("class", "col-lg-auto p-2")
    << run_step_but;

  control_div.Div("button-row")
    << UI::Div("run-col").SetAttr("class", "col-lg-auto p-2")
    << run_toggle_but;

  control_div.Div("button-row")
    << UI::Div("render-frequency-col").SetAttr("class", "col-lg-auto p-2")
    << UI::Div("render-wrapper").SetAttr("class", "input-group input-group-lg")
    << UI::Div("render-input-prepend").SetAttr("class", "input-group-prepend")
    << UI::Div("render-input-prepend-text").SetAttr("class", "input-group-text")
    << "Render every";

  control_div.Div("render-wrapper")
    << UI::Input([this](std::string in) {
                  const size_t val = std::stoul(in);
                  draw_frequency = val;
                 },
                 "number",
                 "",
                 "render-frequency")
        .Checker([this](std::string in) { return CheckInputSize_t(in); })
        .Value((double)draw_frequency)
        .Min(1)
        .Max(std::numeric_limits<size_t>::max())
        .Step(1)
        .SetAttr("class", "form-control")
        .SetCSS("min-width", "96px");

   control_div.Div("render-wrapper")
    << UI::Div("input-group-append").SetAttr("class", "input-group-append")
    << UI::Text().SetAttr("class", "input-group-text")
    << "th generation";

  // ---- Setup config view interface ----
  std::cout << "Setup config interface.."<< std::endl;
  SetupConfigInterface();
  std::cout << "done config interface setup"<< std::endl;

  // ---- Setup world view interface ----

  // Setup world view.
  // - environment -
  world_div << UI::Div("env-canvas-row").SetAttr("class", "row mb-2");
  world_div.Div("env-canvas-row")
    << UI::Div("env-canvas-col").SetAttr("class", "col");

  // - population -
  world_div << UI::Div("pop-canvas-row").SetAttr("class", "row");

  world_div.Div("pop-canvas-row")
    << UI::Div("pop-canvas-col").SetAttr("class", "col")
    << UI::Div("emp-pop-vis");


  // Initial world configuration + pop canvas configuration.
  Setup(); // Call world setup


  // ---- Wire up event handlers ----

  UI::OnDocumentReady([this]() {
    std::cout << "-- OnDocumentReady (open) --"<<std::endl;

    // Configure
    pop_vis.Setup(*this);
    pop_vis.Start();

    emp::OnResize([this]() {
      std::cout << "Resize?" << std::endl;
      RedrawPopulation(false);
      RedrawEnvironment();
    });

    RedrawPopulation(); // Finally, go ahead and draw initial population.
    RedrawEnvironment();
    std::cout << "-- OnDocumentReady (close) --"<<std::endl;
  });


}

void AagosWebInterface::DoFrame() {
  // std::cout << "Frame!" << std::endl;
  RunStep();
  if (GetUpdate() % draw_frequency == 0) {
    RedrawPopulation();
    RedrawEnvironment();
  }
}

void AagosWebInterface::RedrawPopulation(bool update_data) {
  pop_vis.DrawPop(*this, update_data);
}

void AagosWebInterface::RedrawEnvironment() {
  pop_vis.DrawGradientEnv(*this, true);
}

void AagosWebInterface::SetupConfigInterface() {
  config_general_div.Clear();

  // --- Initialize relevant configuration options ---
  // -- general --
  config_input_elements["CHANGE_MAGNITUDE"]={[](std::string){ return true; }, "", ""};
  auto & chg_magnitude_input = config_input_elements["CHANGE_MAGNITUDE"];
  chg_magnitude_input.Type("number");
  chg_magnitude_input.Label("");
  chg_magnitude_input.Min(0);
  chg_magnitude_input.Max(std::numeric_limits<size_t>::max());
  chg_magnitude_input.Step(1);
  chg_magnitude_input.Checker([this](std::string in) { return CheckInputSize_t(in); });
  chg_magnitude_input.Callback([](std::string in) { std::cout << "CHG_MAG" << std::endl; });
  chg_magnitude_input.Value(GetConfig().CHANGE_MAGNITUDE());

  // config_input_elements["CHANGE_FREQUENCY"]={input_callback_fun_t(), "", ""};
  // config_input_elements["POP_SIZE"]={input_callback_fun_t(), "", ""};
  // config_input_elements["MAX_GENS"]={input_callback_fun_t(), "", ""};
  // config_input_elements["SEED"]={input_callback_fun_t(), "", ""};
  // config_input_elements["TOURNAMENT_SIZE"]={input_callback_fun_t(), "", ""};
  // config_input_elements["GRADIENT_MODEL"]={input_callback_fun_t(), "", ""};

  // // // -- genetic architecture --
  // config_input_elements["NUM_BITS_INPUT"]={input_callback_fun_t(), "", ""};
  // config_input_elements["NUM_GENES_INPUT"]={input_callback_fun_t(), "", ""};
  // config_input_elements["GENE_SIZE_INPUT"]={input_callback_fun_t(), "", ""};
  // config_input_elements["MAX_SIZE_INPUT"]={input_callback_fun_t(), "", ""};
  // config_input_elements["MIN_SIZE_INPUT"]={input_callback_fun_t(), "", ""};

  // // // -- Mutations --
  // config_input_elements["GENE_MOVE_PROB"]={input_callback_fun_t(), "", ""};
  // config_input_elements["BIT_FLIP_PROB"]={input_callback_fun_t(), "", ""};
  // config_input_elements["BIT_INS_PROB"]={input_callback_fun_t(), "", ""};
  // config_input_elements["BIT_DEL_PROB"]={input_callback_fun_t(), "", ""};

  // // // -- phase 2 --
  // config_input_elements["PHASE_2_ACTIVE"]={input_callback_fun_t(), "", ""};
  // config_input_elements["PHASE_2_CHANGE_MAGNITUDE"]={input_callback_fun_t(), "", ""};
  // config_input_elements["PHASE_2_CHANGE_FREQUENCY"]={input_callback_fun_t(), "", ""};
  // config_input_elements["PHASE_2_MAX_GENS"]={input_callback_fun_t(), "", ""};
  // config_input_elements["PHASE_2_TOURNAMENT_SIZE"]={input_callback_fun_t(), "", ""};
  // config_input_elements["PHASE_2_GENE_MOVE_PROB"]={input_callback_fun_t(), "", ""};
  // config_input_elements["PHASE_2_BIT_FLIP_PROB"]={input_callback_fun_t(), "", ""};
  // config_input_elements["PHASE_2_BIT_INS_PROB"]={input_callback_fun_t(), "", ""};
  // config_input_elements["PHASE_2_BIT_DEL_PROB"]={input_callback_fun_t(), "", ""};

  // -- Add inputs to interface --

  config_general_div
    << UI::Element("ul", "general-config-ul")
        .SetAttr("class", "list-group list-group-flush")
    << UI::Element("li", "change-magnitude-config-li")
        .SetAttr("class", "list-group-item")
    << config_input_elements["CHANGE_MAGNITUDE"];

}

#endif