#ifndef AAGOS_WEB_D3_H
#define AAGOS_WEB_D3_H

/**
 * TODOs
 * - [ ] add 'phenotype evaluated flag' => use for tooltips instead of current horrible implementation
*/


#include "emp/web/web.hpp"
#include "emp/web/Input.hpp"
#include "emp/web/Element.hpp"
#include "emp/tools/string_utils.hpp"

#include <unordered_map>
#include <sstream>

#include "AagosWorld.h"
#include "AagosConfig.h"
#include "AagosOrg.h"
#include "AagosPopulationVisualization.h"

namespace UI = emp::web;

std::unordered_map<std::string, std::string> config_tooltips = {
  {"CHANGE_MAGNITUDE", "When the environment changes, how many bits are flipped across all gene targets? (0 for no change)"},
  {"CHANGE_FREQUENCY", "How many generations elapse between environmental changes? (0 for no change)"},
  {"POP_SIZE", "How many digital genomes are in the population?"},
  {"MAX_GENS", "How many generations should we evolve the population? Generations are synchronous. All genomes are evaluated and selected for asexual reproduction."},
  {"TOURNAMENT_SIZE", "We use tournament selection to determine which genomes reproduce. Tournament size controls the strength of selection (higher is stronger, lower is weaker; 1 results in pure genetic drift)"},
  {"SEED", "Random number generator seed."},

  {"NUM_GENES", "How many genes does each genome have?"},
  {"GENE_SIZE", "Gene size (in bits)"},
  {"MIN_SIZE", "Minimum genome size (in bits)"},
  {"MAX_SIZE", "Maximum genome size (in bits)"},
  {"NUM_BITS", "Starting genome size (in bits)"},

  {"GENE_MOVE_PROB", "Per-gene probability of mutating gene location"},
  {"BIT_FLIP_PROB", "Per-bit probability of bit substitutions"},
  {"BIT_INS_PROB", "Per-bit probability of a single bit insertion"},
  {"BIT_DEL_PROB", "Per-bit probability of a single bit deletion"},
  {"APPLY_BIT_MUTS_PER_GENE", "Apply bit mutations to coding sites at a per-gene-per-site rate instead of just per-site (overlapping genes no longer reduce mutational load on coding region)."},

  {"PHASE_2_ACTIVE", "After MAX_GENS number of generations have elapsed, should run continue to a second phase with new parameters (below)?"},
  {"PHASE_2_CHANGE_MAGNITUDE", "Phase 2 - When the environment changes, how many bits are flipped across all gene targets? (0 for no change)"},
  {"PHASE_2_CHANGE_FREQUENCY", "Phase 2 - How many generations elapse between environmental changes? (0 for no change)"},
  {"PHASE_2_MAX_GENS", "Phase 2 - How many generations should we evolve the population? Generations are synchronous. All genomes are evaluated and selected for asexual reproduction."},
  {"PHASE_2_TOURNAMENT_SIZE", "Phase 2 - We use tournament selection to determine which genomes reproduce. Tournament size controls the strength of selection (higher is stronger, lower is weaker; 1 results in pure genetic drift)"},
  {"PHASE_2_GENE_MOVE_PROB", "Per-gene probability of mutating gene location"},
  {"PHASE_2_BIT_FLIP_PROB", "Per-bit probability of bit substitutions"},
  {"PHASE_2_BIT_INS_PROB", "Per-bit probability of a single bit insertion"},
  {"PHASE_2_BIT_DEL_PROB", "Per-bit probability of a single bit deletion"}
};

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

// https://stackoverflow.com/questions/29169153/how-do-i-verify-a-string-is-valid-double-even-if-it-has-a-point-in-it
bool is_numeric(const std::string & str) {
    double result = double();
    auto i = std::istringstream(str);
    i >> result;
    return !i.fail() && i.eof();
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
  UI::Document stats_div;
  UI::Document pop_vis_div;
  UI::Document control_div;
  UI::Document config_general_div;
  UI::Document config_genetic_arch_div;
  UI::Document config_mutation_div;
  UI::Document config_phase_2_evo_div;
  UI::Document confirm_exp_config_div;

  UI::Button run_toggle_but;
  UI::Button run_step_but;
  UI::Button run_reset_but;
  UI::Button config_exp_but;
  UI::Button config_apply_but;
  UI::Button confirm_config_exp_but;
  UI::Button draw_mode_toggle_but;

  std::unordered_map<std::string, UI::Input> config_input_elements;

  size_t draw_frequency=64;
  size_t generation_last_drawn=0;
  bool config_mode=false;

  AagosPopulationVisualization pop_vis;

  void RedrawPopulation(bool update_data=true);
  void RedrawEnvironment();
  void ReconfigureWorld();

  void SetupConfigInterface();
  void SetupStatsViewInterface();

  UI::Input & AddConfigInput(UI::Document & div,
                             const std::string & append_to_id,
                             const std::string & config_name,
                             const std::string & type,
                             const input_checker_fun_t & input_checker,
                             const std::string & min_val="",
                             const std::string & max_val="",
                             const std::string & step_val="",
                             const std::string & init_val="",
                             const std::string & config_tooltip="")
  {
    // If config element doesn't exist yet, create a stub.
    if (!emp::Has(config_input_elements, config_name)) {
      config_input_elements[config_name]={[](std::string){ return true; }, "", ""};
    }
    // Configure config input.
    auto & input_elem = config_input_elements[config_name];
    input_elem.Type(type);
    input_elem.Label("");
    input_elem.Checker(input_checker);
    input_elem.Min(min_val);
    input_elem.Max(max_val);
    input_elem.Step(step_val);
    input_elem.Callback([config_name, this](std::string in) {
      std::cout << "Callback: " << config_name << std::endl;
      std::cout << "  Value = " << in << std::endl;
      // pretty sure that the callback only gets triggered if check is passed
      config_input_elements[config_name].Value(in);
    });
    input_elem.Value(init_val);
    input_elem.SetCSS("min-width", "64px");       // todo - will this work for all?
    input_elem.SetAttr("class", "form-control");

    div.Find(append_to_id)
      << UI::Div(config_name + "-config-wrapper").SetAttr("class", "input-group")
      << UI::Div(config_name + "-config-input-prepend").SetAttr("class", "input-group-prepend")
          .SetAttr("data-toggle", "tooltip")
          .SetAttr("data-placement", "top")
          .SetAttr("title", config_tooltip)
      << UI::Div(config_name + "-config-input-prepend-text").SetAttr("class", "input-group-text")
    << config_name;

    if (type == "checkbox") {
      // Do slightly different HTML structure for checkboxes (slightly improves bootstrap styling)
      div.Find(config_name + "-config-wrapper")
        << UI::Div(config_name + "-config-input-append").SetAttr("class", "input-group-append")
        << UI::Div(config_name + "-config-input-append-text").SetAttr("class", "input-group-text")
        << config_input_elements[config_name];
    } else {
      div.Find(config_name + "-config-wrapper")
        << config_input_elements[config_name];
    }

    return config_input_elements[config_name];
  }

  void DisableConfigInputs(bool in_dis=true) {
    for (auto & cfg : config_input_elements) {
      cfg.second.Disabled(in_dis);
    }
  }

  bool CheckInputSize_t(const std::string & in) {
    return (emp::is_digits(in) && in.size() && std::stoi(in) >= 0); // todo - FromString?
  }

  bool CheckInputDouble(const std::string & in) {
    return (in.size() && is_numeric(in));
  }

  bool CheckInputProbability(const std::string & in) {
    if (!CheckInputDouble(in)) return false;
    const double prob = emp::from_string<double>(in);
    return prob >= 0.0 && prob <= 1.0;
  }

  bool CheckInput01(const std::string & in) {
    if (!CheckInputSize_t(in)) return false;
    const size_t val = emp::from_string<size_t>(in);
    return val == 0 || val == 1;
  }

  bool CheckInput_GENE_SIZE(const std::string & in) {
    if (!CheckInputSize_t(in)) return false;
    const size_t val = emp::from_string<size_t>(in);
    const size_t max_size = emp::Has(config_input_elements, "MAX_SIZE")
                                  ? emp::from_string<size_t>(config_input_elements["MAX_SIZE"].GetValue())
                                  : val;
    // const size_t min_size = emp::Has(config_input_elements, "MIN_SIZE")
    //                               ? emp::from_string<size_t>(config_input_elements["MIN_SIZE"].GetValue())
    //                               : val;
    return val <= max_size;
  }

  // todo - make such that this doesn't depend on order of input initialization!
  bool CheckInput_MIN_SIZE(const std::string & in) {
    if (!CheckInputSize_t(in)) return false;
    const size_t val = emp::from_string<size_t>(in);
    const size_t cur_gene_size = emp::from_string<size_t>(config_input_elements["GENE_SIZE"].GetValue());
    const size_t cur_max_size = emp::Has(config_input_elements, "MAX_SIZE")
                                  ? emp::from_string<size_t>(config_input_elements["MAX_SIZE"].GetValue())
                                  : std::numeric_limits<size_t>::max();
    const size_t num_bits = emp::Has(config_input_elements, "NUM_BITS")
                                  ? emp::from_string<size_t>(config_input_elements["NUM_BITS"].GetValue())
                                  : val;
    return val >= cur_gene_size && val <= cur_max_size && val <= num_bits;
  }

  // todo - make such that this doesn't depend on order of input initialization!
  bool CheckInput_MAX_SIZE(const std::string & in) {
    if (!CheckInputSize_t(in)) { return false; }
    const size_t in_num = emp::from_string<size_t>(in);
    const size_t gene_size = emp::from_string<size_t>(config_input_elements["GENE_SIZE"].GetValue());
    const size_t min_size = emp::from_string<size_t>(config_input_elements["MIN_SIZE"].GetValue());
    const size_t num_bits = emp::Has(config_input_elements, "NUM_BITS")
                                  ? emp::from_string<size_t>(config_input_elements["NUM_BITS"].GetValue())
                                  : in_num;
    return (in_num >= gene_size) && (in_num >= min_size) && (in_num >= num_bits);
  }

  // todo - make such that this doesn't depend on order of input initialization!
  bool CheckInput_NUM_BITS(const std::string & in) {
    if (!CheckInputSize_t(in)) return false;
    const size_t in_num = emp::from_string<size_t>(in);
    const size_t gene_size = emp::from_string<size_t>(config_input_elements["GENE_SIZE"].GetValue());
    const size_t min_size = emp::from_string<size_t>(config_input_elements["MIN_SIZE"].GetValue());
    const size_t max_size = emp::from_string<size_t>(config_input_elements["MAX_SIZE"].GetValue());
    return (in_num >= gene_size) && (in_num <= max_size) && (in_num >= min_size);
  }

public:
  AagosWebInterface(config_t & cfg)
    : AagosWorld(cfg),
      world_div("emp_world_view"),
      stats_div("emp_world_stats_view"),
      pop_vis_div("emp_pop_vis_view"),
      control_div("emp_controls_view"),
      config_general_div("emp_view_config_general"),
      config_genetic_arch_div("emp_view_config_genetic_architecture"),
      config_mutation_div("emp_view_config_mutation"),
      config_phase_2_evo_div("emp_view_config_phase_2_evolution"),
      confirm_exp_config_div("emp_config_exp_confirmation"),
      pop_vis("emp_pop_vis_view")
      // pop_vis("emp-pop-vis")
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

  // emp::JSWrap([this](){ ReconfigureWorld(); }, "reconfigure_world");
  emp::JSWrap([this]() {

    ReconfigureWorld();

    config_exp_but.SetAttr("class", "btn btn-block btn-lg btn-danger");
    config_apply_but.SetAttr("class", "d-none");
    config_apply_but.SetLabel("Apply configuration");

    run_toggle_but.SetDisabled(config_mode);
    run_step_but.SetDisabled(config_mode);
    run_reset_but.SetDisabled(config_mode);

    DisableConfigInputs(!config_mode);
  }, "do_apply_config"); // Yikes, yucking up empirical's 'global' js namespace

  emp::JSWrap([this]() {
    ReconfigureWorld();
    run_reset_but.SetDisabled(active);
    run_toggle_but.SetDisabled(active);
    run_step_but.SetDisabled(active);
    config_exp_but.SetDisabled(active);
  }, "do_reset_world");

  // ---- Setup interface control buttons ----
  // Run button setup.
  // Manually create run/pause button to have full control over button's callback.
  // NOTE - here's an example where Empirical's web tools are overly restrictive.
  run_toggle_but = UI::Button([this]() {
    ToggleActive();
    run_toggle_but.SetLabel(active ? "Stop" : "Run"); // TODO - use icons
    run_reset_but.SetDisabled(active);
    run_step_but.SetDisabled(active);
    config_exp_but.SetDisabled(active);
  }, "Run", "run-toggle-button");
  run_toggle_but.SetAttr("class", "btn btn-block btn-lg btn-primary");

  // Step button setup.
  run_step_but = GetStepButton("run-step-button");
  run_step_but.SetAttr("class", "btn btn-block btn-lg btn-primary");

  // Reset run button setup.
  run_reset_but = UI::Button([]() {
    EM_ASM({
      // Hi-jack button's html to show dope spinner while doing World's Setup (which can take a min)
      $("#run-reset-button").html('<div class="d-flex align-items-center justify-content-center"><span class="spinner-border spinner-border-sm mr-1" role="status"></span>Restarting...</div>');
      setTimeout(function() {
        emp.do_reset_world();
      });
    }, 42);
  }, "Restart", "run-reset-button");
  run_reset_but.SetAttr("class", "btn btn-block btn-lg btn-warning");

  // Config experiment button setup.
  config_exp_but = UI::Button([]() { ; }, "Configure Experiment", "config-exp-button");
  config_exp_but.SetAttr("class", "btn btn-block btn-lg btn-danger");
  config_exp_but.SetAttr("data-toggle", "modal");
  config_exp_but.SetAttr("data-target", "#config-exp-modal");

  config_apply_but = UI::Button([this]() {
    // In config mode, so we need to apply configuration.
    config_mode = false;

    EM_ASM({
      // Hi-jack config-apply-button's html to show dope spinner while doing World's Setup (which can take a min)
      $("#config-apply-button").html('<div class="d-flex align-items-center justify-content-center"><span class="spinner-border spinner-border-sm mr-1" role="status"></span>Applying configuration...</div>');
      setTimeout(function() {
        emp.do_apply_config();
        // $("#config-loading-spinner").attr("class", "spinner-border d-none");
      });
    }, 42);

  }, "Apply configuration", "config-apply-button");
  config_apply_but.SetAttr("class", "d-none");


  confirm_config_exp_but = UI::Button([this]() {
    config_mode = true;

    config_exp_but.SetAttr("class", "d-none");
    config_apply_but.SetAttr("class", "btn btn-block btn-lg btn-danger");

    run_toggle_but.SetDisabled(config_mode);
    run_step_but.SetDisabled(config_mode);
    run_reset_but.SetDisabled(config_mode);
    DisableConfigInputs(!config_mode);

    // TODO - clear d3 vis
    pop_vis.Clear();

  }, "Confirm", "confirm-config-exp-button");
  confirm_config_exp_but.SetAttr("class", "btn btn-danger");
  confirm_config_exp_but.SetAttr("data-dismiss", "modal");
  confirm_exp_config_div << confirm_config_exp_but;

  draw_mode_toggle_but = UI::Button([this]() {
    if (pop_vis.IsDrawModeFullPop()) {
      pop_vis.SetDrawModeMaxFit();
      draw_mode_toggle_but.SetLabel("Draw Full Population");
    } else if (pop_vis.IsDrawModeMaxFit()) {
      pop_vis.SetDrawModeFullPop();
      draw_mode_toggle_but.SetLabel("Draw Max Fitness Organism");
    }
    if (!config_mode) {
      RedrawPopulation(true); // update data, disable tooltips
      RedrawEnvironment();
    }
  }, "Draw Full Population", "population-draw-mode-toggle-button");
  draw_mode_toggle_but.SetAttr("class", "btn btn-block btn-lg btn-secondary");

  // Add buttons to controls view.
  control_div << UI::Div("button-row").SetAttr("class", "row justify-content-md-center");

  control_div.Div("button-row")
    << UI::Div("step-col").SetAttr("class", "col-lg-auto p-2")
    << run_step_but;

  control_div.Div("button-row")
    << UI::Div("run-col").SetAttr("class", "col-lg-auto p-2")
    << run_toggle_but;

  control_div.Div("button-row")
    << UI::Div("reset-col").SetAttr("class", "col-lg-auto p-2")
    << run_reset_but;

  control_div.Div("button-row")
    << UI::Div("config-exp-col").SetAttr("class", "col-lg-auto p-2")
    << config_exp_but;
  control_div.Div("config-exp-col")
    << config_apply_but;

  control_div.Div("button-row")
    << UI::Div("draw-mode-col").SetAttr("class", "col-lg-auto p-2")
    << draw_mode_toggle_but;

  control_div.Div("button-row")
    << UI::Div("render-frequency-col").SetAttr("class", "col-lg-auto p-2")
    // << UI::Div("render-wrapper").SetAttr("class", "input-group input-group-lg")
    << UI::Div("render-wrapper")
        .SetAttr("class", "input-group input-group-md")
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
        .SetCSS("min-width", "64px");

  control_div.Div("render-wrapper")
    << UI::Div("input-group-append").SetAttr("class", "input-group-append")
    << UI::Text().SetAttr("class", "input-group-text")
    << "generations";


  // ---- Setup config view interface ----
  std::cout << "Setup config interface.."<< std::endl;
  SetupConfigInterface();
  std::cout << "done config interface setup"<< std::endl;

  // ---- Setup world view interface ----

  // Setup world view.
  // - environment -
  world_div << UI::Div("hud-row").SetAttr("class", "row justify-content-center");
  world_div.Div("hud-row")
    << UI::Div("generation-counter-col").SetAttr("class", "col-sm-auto p-1")
    << UI::Element("h4", "")
    << UI::Element("span", "generation-counter-badge").SetAttr("class", "badge badge-secondary")
    << "Generation: ";
  world_div.Div("generation-counter-badge")
    << UI::Live([this]() { return GetUpdate(); });

  world_div.Div("hud-row")
    << UI::Div("phase-counter-col").SetAttr("class", "col-sm-auto p-1")
    << UI::Element("h4", "")
    << UI::Element("span", "phase-counter-badge").SetAttr("class", "badge badge-secondary")
    << "Experiment phase: ";
  world_div.Div("phase-counter-badge")
    << UI::Live([this]() { return cur_phase; });

  world_div.Div("hud-row")
    << UI::Div("last-rendered-col").SetAttr("class", "col-sm-auto p-1")
    << UI::Element("h4", "")
    << UI::Element("span", "generation-last-rendered-counter-badge").SetAttr("class", "badge badge-secondary")
    << "Generation last drawn: ";
  world_div.Div("generation-last-rendered-counter-badge")
    << UI::Live([this]() { return generation_last_drawn; });

  world_div << UI::Element("hr").SetAttr("class", "mt-1");

  // --- Configure stats ---
  SetupStatsViewInterface();

  // Initial world configuration + pop canvas configuration.
  Setup(); // Call world setup

  // ---- Wire up event handlers ----
  UI::OnDocumentReady([this]() {
    std::cout << "-- OnDocumentReady (open) --"<<std::endl;

    // Configure
    pop_vis.Setup(*this);

    emp::OnResize([this]() {
      std::cout << "Resize?" << std::endl;
      RedrawPopulation(false);
      RedrawEnvironment();
    });

    // RedrawPopulation(true, false); // Finally, go ahead and draw initial population.
    // RedrawEnvironment();
    DoFrame();

    // Enable tooltips!
    EM_ASM({
      $(function () {
        $('[data-toggle="tooltip"]').tooltip(({container: 'body', boundary: 'window'}));
      });
    });

    std::cout << "-- OnDocumentReady (close) --"<<std::endl;
  });

}

void AagosWebInterface::DoFrame() {
  // std::cout << "Frame!" << std::endl;
  RunStep(false);
  if (   (GetUpdate() % draw_frequency == 0)
      || (GetUpdate() == config.MAX_GENS())
      || (GetUpdate() == TOTAL_GENS))
  {
    RedrawPopulation(true);
    RedrawEnvironment();
  }

  // Do all the checks, etc
  if ((cur_phase==0) && (GetUpdate() >= config.MAX_GENS())) {
    // Trigger phase 1 || stop!
    if (config.PHASE_2_ACTIVE()) {
      ActivateEvoPhaseTwo();
    } else {
      // Stop!
      ToggleActive();
      run_toggle_but.SetLabel("Run");
      run_toggle_but.SetDisabled(true);
      run_step_but.SetDisabled(true);
      config_exp_but.SetDisabled(false);
      run_reset_but.SetDisabled(false);
    }

  } else if (cur_phase == 1 && GetUpdate() >= TOTAL_GENS) {
    // Trigger stop!
    ToggleActive();
    run_toggle_but.SetLabel("Run");
    run_toggle_but.SetDisabled(true);
    run_step_but.SetDisabled(true);
    config_exp_but.SetDisabled(false);
    run_reset_but.SetDisabled(false);
  }
  world_div.Redraw();
  stats_div.Redraw();

  // std::cout << "Most fit ID = " << GetMostFitID() << std::endl;
  // std::cout << " - Fitness = " << GetOrg(GetMostFitID()).GetPhenotype().fitness << std::endl;
  // std::cout << " - Gene contributions = ";
  // emp::Print(GetOrg(GetMostFitID()).GetPhenotype().gene_fitness_contributions);
  // std::cout << std::endl;
  // std::cout << " - Bits = " << GetOrg(GetMostFitID()).GetGenome().bits << std::endl;
  // std::cout << " - Gene starts = ";
  // emp::Print(GetOrg(GetMostFitID()).GetGenome().gene_starts);
  // std::cout << std::endl;

  // Manually do the world update + clear world cache.
  AdvanceWorld();
}

void AagosWebInterface::RedrawPopulation(bool update_data/*=true*/) {
  if (update_data) {
    generation_last_drawn = GetUpdate();
  }
  pop_vis.DrawPop(*this, update_data);
}

void AagosWebInterface::RedrawEnvironment() {
  pop_vis.DrawGradientEnv(*this, true);
}

void AagosWebInterface::ReconfigureWorld() {
  std::cout << "--- reconfigure world ---" << std::endl;

  // Loop over config inputs, reconfiguring.
  for (auto & cfg : config_input_elements) {
    config.Set(cfg.first, cfg.second.GetValue());
  }
  // Setup the world again...
  Setup();
  pop_vis.Setup(*this); // Re-configure pop_vis
  GetFitnessDataNode()->Reset();

  DoFrame();
  // RedrawPopulation(true, false); // Finally, go ahead and draw initial population.
  // RedrawEnvironment();
  // world_div.Redraw();
  // stats_div.Redraw();
  std::cout << "--- done reconfiguring world ---" << std::endl;
}

void AagosWebInterface::SetupConfigInterface() {
  config_general_div.Clear();

  // --- Initialize relevant configuration options ---
  // -- general --

  // --- General Configuration Settings ---
  config_general_div
    << UI::Element("ul", "general-config-ul")
        .SetAttr("class", "list-group list-group-flush");
        // .SetCSS("overflow-x","scroll");

  // CHANGE_MAGNITUDE
  config_general_div.Find("general-config-ul")
    << UI::Element("li", "CHANGE_MAGNITUDE-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_general_div,
                /* append_to_id   =*/ "CHANGE_MAGNITUDE-config-li",
                /* config_name    =*/ "CHANGE_MAGNITUDE",
                /* type           =*/ "number",
                /* checker        =*/ [this](std::string in) { return CheckInputSize_t(in); },
                /* min_val        =*/ "0",
                /* max_val        =*/ emp::to_string(std::numeric_limits<size_t>::max()),
                /* step_val       =*/ "1",
                /* init_val       =*/ emp::to_string(GetConfig().CHANGE_MAGNITUDE()),
                // /* config_tooltip =*/ GetConfig()["CHANGE_MAGNITUDE"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["CHANGE_MAGNITUDE"]);

  // CHANGE_FREQUENCY
  config_general_div.Find("general-config-ul")
    << UI::Element("li", "CHANGE_FREQUENCY-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_general_div,
                /* append_to_id =*/ "CHANGE_FREQUENCY-config-li",
                /* config_name  =*/ "CHANGE_FREQUENCY",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInputSize_t(in); },
                /* min_val      =*/ "0",
                /* max_val      =*/ emp::to_string(std::numeric_limits<size_t>::max()),
                /* step_val     =*/ "1",
                /* init_val     =*/ emp::to_string(GetConfig().CHANGE_FREQUENCY()),
                // /* config_tooltip =*/ GetConfig()["CHANGE_FREQUENCY"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["CHANGE_FREQUENCY"]);

  // POP_SIZE
  config_general_div.Find("general-config-ul")
    << UI::Element("li", "POP_SIZE-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_general_div,
                /* append_to_id =*/ "POP_SIZE-config-li",
                /* config_name  =*/ "POP_SIZE",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInputSize_t(in); },
                /* min_val      =*/ "1",
                /* max_val      =*/ emp::to_string(std::numeric_limits<size_t>::max()),
                /* step_val     =*/ "1",
                /* init_val     =*/ emp::to_string(GetConfig().POP_SIZE()),
                // /* config_tooltip =*/ GetConfig()["POP_SIZE"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["POP_SIZE"]);


  // MAX_GENS
  config_general_div.Find("general-config-ul")
    << UI::Element("li", "MAX_GENS-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_general_div,
                /* append_to_id =*/ "MAX_GENS-config-li",
                /* config_name  =*/ "MAX_GENS",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInputSize_t(in); },
                /* min_val      =*/ "0",
                /* max_val      =*/ emp::to_string(std::numeric_limits<size_t>::max()),
                /* step_val     =*/ "1",
                /* init_val     =*/ emp::to_string(GetConfig().MAX_GENS()),
                // /* config_tooltip =*/ GetConfig()["MAX_GENS"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["MAX_GENS"]);

  // TOURNAMENT_SIZE
  config_general_div.Find("general-config-ul")
    << UI::Element("li", "TOURNAMENT_SIZE-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_general_div,
                /* append_to_id =*/ "TOURNAMENT_SIZE-config-li",
                /* config_name  =*/ "TOURNAMENT_SIZE",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInputSize_t(in); },
                /* min_val      =*/ "1",
                /* max_val      =*/ emp::to_string(std::numeric_limits<size_t>::max()),
                /* step_val     =*/ "1",
                /* init_val     =*/ emp::to_string(GetConfig().TOURNAMENT_SIZE()),
                // /* config_tooltip =*/ GetConfig()["TOURNAMENT_SIZE"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["TOURNAMENT_SIZE"]);

  // SEED
  config_general_div.Find("general-config-ul")
    << UI::Element("li", "SEED-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_general_div,
                /* append_to_id =*/ "SEED-config-li",
                /* config_name  =*/ "SEED",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInputSize_t(in); },
                /* min_val      =*/ "0",
                /* max_val      =*/ emp::to_string(std::numeric_limits<size_t>::max()),
                /* step_val     =*/ "1",
                /* init_val     =*/ emp::to_string(GetConfig().SEED()),
                // /* config_tooltip =*/ GetConfig()["SEED"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["SEED"]);


  // -- genetic architecture --
  config_genetic_arch_div
    << UI::Element("ul", "genetic-arch-config-ul")
        .SetAttr("class", "list-group list-group-flush")
        .SetCSS("overflow-x","scroll");

  // NUM_GENES
  config_genetic_arch_div.Find("genetic-arch-config-ul")
    << UI::Element("li", "NUM_GENES-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_genetic_arch_div,
                /* append_to_id =*/ "NUM_GENES-config-li",
                /* config_name  =*/ "NUM_GENES",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInputSize_t(in); },
                /* min_val      =*/ "1",
                /* max_val      =*/ emp::to_string(std::numeric_limits<size_t>::max()),
                /* step_val     =*/ "1",
                /* init_val     =*/ emp::to_string(GetConfig().NUM_GENES()),
                // /* config_tooltip =*/ GetConfig()["NUM_GENES"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["NUM_GENES"]);

  // GENE_SIZE
  config_genetic_arch_div.Find("genetic-arch-config-ul")
    << UI::Element("li", "GENE_SIZE-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_genetic_arch_div,
                /* append_to_id =*/ "GENE_SIZE-config-li",
                /* config_name  =*/ "GENE_SIZE",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInput_GENE_SIZE(in); },
                /* min_val      =*/ "1",
                /* max_val      =*/ emp::to_string(std::numeric_limits<size_t>::max()),
                /* step_val     =*/ "1",
                /* init_val     =*/ emp::to_string(GetConfig().GENE_SIZE()),
                // /* config_tooltip =*/ GetConfig()["GENE_SIZE"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["GENE_SIZE"]);

  // MIN_SIZE
  config_genetic_arch_div.Find("genetic-arch-config-ul")
    << UI::Element("li", "MIN_SIZE-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_genetic_arch_div,
                /* append_to_id =*/ "MIN_SIZE-config-li",
                /* config_name  =*/ "MIN_SIZE",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInput_MIN_SIZE(in); },
                /* min_val      =*/ "1",
                /* max_val      =*/ emp::to_string(std::numeric_limits<size_t>::max()),
                /* step_val     =*/ "1",
                /* init_val     =*/ emp::to_string(GetConfig().MIN_SIZE()),
                // /* config_tooltip =*/ GetConfig()["MIN_SIZE"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["MIN_SIZE"]);

  // MAX_SIZE
  config_genetic_arch_div.Find("genetic-arch-config-ul")
    << UI::Element("li", "MAX_SIZE-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_genetic_arch_div,
                /* append_to_id =*/ "MAX_SIZE-config-li",
                /* config_name  =*/ "MAX_SIZE",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInput_MAX_SIZE(in); },
                /* min_val      =*/ "1",
                /* max_val      =*/ emp::to_string(std::numeric_limits<size_t>::max()),
                /* step_val     =*/ "1",
                /* init_val     =*/ emp::to_string(GetConfig().MAX_SIZE()),
                // /* config_tooltip =*/ GetConfig()["MAX_SIZE"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["MAX_SIZE"]);

  // NUM_BITS
  config_genetic_arch_div.Find("genetic-arch-config-ul")
    << UI::Element("li", "NUM_BITS-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_genetic_arch_div,
                /* append_to_id =*/ "NUM_BITS-config-li",
                /* config_name  =*/ "NUM_BITS",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInput_NUM_BITS(in); },
                /* min_val      =*/ "0",
                /* max_val      =*/ emp::to_string(std::numeric_limits<size_t>::max()),
                /* step_val     =*/ "1",
                /* init_val     =*/ emp::to_string(GetConfig().NUM_BITS()),
                // /* config_tooltip =*/ GetConfig()["NUM_BITS"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["NUM_BITS"]);

  // -- Mutations --
  config_mutation_div
    << UI::Element("ul", "mutation-config-ul")
        .SetAttr("class", "list-group list-group-flush")
        .SetCSS("overflow-x","scroll");

  // GENE_MOVE_PROB
  config_mutation_div.Find("mutation-config-ul")
    << UI::Element("li", "GENE_MOVE_PROB-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_mutation_div,
                /* append_to_id =*/ "GENE_MOVE_PROB-config-li",
                /* config_name  =*/ "GENE_MOVE_PROB",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInputProbability(in); },
                /* min_val      =*/ "0",
                /* max_val      =*/ "1.0",
                /* step_val     =*/ "any",
                /* init_val     =*/ emp::to_string(GetConfig().GENE_MOVE_PROB()),
                // /* config_tooltip =*/ GetConfig()["GENE_MOVE_PROB"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["GENE_MOVE_PROB"]);

  // BIT_FLIP_PROB
  config_mutation_div.Find("mutation-config-ul")
    << UI::Element("li", "BIT_FLIP_PROB-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_mutation_div,
                /* append_to_id =*/ "BIT_FLIP_PROB-config-li",
                /* config_name  =*/ "BIT_FLIP_PROB",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInputProbability(in); },
                /* min_val      =*/ "0",
                /* max_val      =*/ "1.0",
                /* step_val     =*/ "any",
                /* init_val     =*/ emp::to_string(GetConfig().BIT_FLIP_PROB()),
                // /* config_tooltip =*/ GetConfig()["BIT_FLIP_PROB"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["BIT_FLIP_PROB"]);

  // BIT_INS_PROB
  config_mutation_div.Find("mutation-config-ul")
    << UI::Element("li", "BIT_INS_PROB-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_mutation_div,
                /* append_to_id =*/ "BIT_INS_PROB-config-li",
                /* config_name  =*/ "BIT_INS_PROB",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInputProbability(in); },
                /* min_val      =*/ "0",
                /* max_val      =*/ "1.0",
                /* step_val     =*/ "any",
                /* init_val     =*/ emp::to_string(GetConfig().BIT_INS_PROB()),
                // /* config_tooltip =*/ GetConfig()["BIT_INS_PROB"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["BIT_INS_PROB"]);

  // BIT_DEL_PROB
  config_mutation_div.Find("mutation-config-ul")
    << UI::Element("li", "BIT_DEL_PROB-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_mutation_div,
                /* append_to_id =*/ "BIT_DEL_PROB-config-li",
                /* config_name  =*/ "BIT_DEL_PROB",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInputProbability(in); },
                /* min_val      =*/ "0",
                /* max_val      =*/ "1.0",
                /* step_val     =*/ "any",
                /* init_val     =*/ emp::to_string(GetConfig().BIT_DEL_PROB()),
                // /* config_tooltip =*/ GetConfig()["BIT_DEL_PROB"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["BIT_DEL_PROB"]);

  // APPLY_BIT_MUTS_PER_GENE
  config_mutation_div.Find("mutation-config-ul")
    << UI::Element("li", "APPLY_BIT_MUTS_PER_GENE-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
  AddConfigInput(config_mutation_div,
                /* append_to_id =*/ "APPLY_BIT_MUTS_PER_GENE-config-li",
                /* config_name  =*/ "APPLY_BIT_MUTS_PER_GENE",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInput01(in); },
                /* min_val      =*/ "0",
                /* max_val      =*/ "1",
                /* step_val     =*/ "1",
                /* init_val     =*/ emp::to_string(GetConfig().APPLY_BIT_MUTS_PER_GENE()),
                /* config_tooltip =*/ config_tooltips["APPLY_BIT_MUTS_PER_GENE"]);

  // -- phase 2 --
  config_phase_2_evo_div
    << UI::Element("ul", "phase-2-config-ul")
        .SetAttr("class", "list-group list-group-flush")
        .SetCSS("overflow-x","scroll");

  // PHASE_2_ACTIVE
  config_phase_2_evo_div.Find("phase-2-config-ul")
    << UI::Element("li", "PHASE_2_ACTIVE-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_phase_2_evo_div,
                /* append_to_id =*/ "PHASE_2_ACTIVE-config-li",
                /* config_name  =*/ "PHASE_2_ACTIVE",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInput01(in); },
                /* min_val      =*/ "0",
                /* max_val      =*/ "1",
                /* step_val     =*/ "1",
                /* init_val     =*/ emp::to_string(GetConfig().PHASE_2_ACTIVE()),
                // /* config_tooltip =*/ GetConfig()["PHASE_2_ACTIVE"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["PHASE_2_ACTIVE"]);

  // PHASE_2_CHANGE_MAGNITUDE
  config_phase_2_evo_div.Find("phase-2-config-ul")
    << UI::Element("li", "PHASE_2_CHANGE_MAGNITUDE-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_phase_2_evo_div,
                /* append_to_id =*/ "PHASE_2_CHANGE_MAGNITUDE-config-li",
                /* config_name  =*/ "PHASE_2_CHANGE_MAGNITUDE",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInputSize_t(in); },
                /* min_val      =*/ "0",
                /* max_val      =*/ emp::to_string(std::numeric_limits<size_t>::max()),
                /* step_val     =*/ "1",
                /* init_val     =*/ emp::to_string(GetConfig().PHASE_2_CHANGE_MAGNITUDE()),
                // /* config_tooltip =*/ GetConfig()["PHASE_2_CHANGE_MAGNITUDE"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["PHASE_2_CHANGE_MAGNITUDE"]);

  // PHASE_2_CHANGE_FREQUENCY
  config_phase_2_evo_div.Find("phase-2-config-ul")
    << UI::Element("li", "PHASE_2_CHANGE_FREQUENCY-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_phase_2_evo_div,
                /* append_to_id =*/ "PHASE_2_CHANGE_FREQUENCY-config-li",
                /* config_name  =*/ "PHASE_2_CHANGE_FREQUENCY",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInputSize_t(in); },
                /* min_val      =*/ "0",
                /* max_val      =*/ emp::to_string(std::numeric_limits<size_t>::max()),
                /* step_val     =*/ "1",
                /* init_val     =*/ emp::to_string(GetConfig().PHASE_2_CHANGE_FREQUENCY()),
                // /* config_tooltip =*/ GetConfig()["PHASE_2_CHANGE_FREQUENCY"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["PHASE_2_CHANGE_FREQUENCY"]);

  // PHASE_2_MAX_GENS
  config_phase_2_evo_div.Find("phase-2-config-ul")
    << UI::Element("li", "PHASE_2_MAX_GENS-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_phase_2_evo_div,
                /* append_to_id =*/ "PHASE_2_MAX_GENS-config-li",
                /* config_name  =*/ "PHASE_2_MAX_GENS",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInputSize_t(in); },
                /* min_val      =*/ "0",
                /* max_val      =*/ emp::to_string(std::numeric_limits<size_t>::max()),
                /* step_val     =*/ "1",
                /* init_val     =*/ emp::to_string(GetConfig().PHASE_2_MAX_GENS()),
                // /* config_tooltip =*/ GetConfig()["PHASE_2_MAX_GENS"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["PHASE_2_MAX_GENS"]);

  // PHASE_2_TOURNAMENT_SIZE
  config_phase_2_evo_div.Find("phase-2-config-ul")
    << UI::Element("li", "PHASE_2_TOURNAMENT_SIZE-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_phase_2_evo_div,
                /* append_to_id =*/ "PHASE_2_TOURNAMENT_SIZE-config-li",
                /* config_name  =*/ "PHASE_2_TOURNAMENT_SIZE",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInputSize_t(in); },
                /* min_val      =*/ "1",
                /* max_val      =*/ emp::to_string(std::numeric_limits<size_t>::max()),
                /* step_val     =*/ "1",
                /* init_val     =*/ emp::to_string(GetConfig().PHASE_2_TOURNAMENT_SIZE()),
                // /* config_tooltip =*/ GetConfig()["PHASE_2_TOURNAMENT_SIZE"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["PHASE_2_TOURNAMENT_SIZE"]);

  // PHASE_2_GENE_MOVE_PROB
  config_phase_2_evo_div.Find("phase-2-config-ul")
    << UI::Element("li", "PHASE_2_GENE_MOVE_PROB-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_phase_2_evo_div,
                /* append_to_id =*/ "PHASE_2_GENE_MOVE_PROB-config-li",
                /* config_name  =*/ "PHASE_2_GENE_MOVE_PROB",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInputProbability(in); },
                /* min_val      =*/ "0",
                /* max_val      =*/ "1.0",
                /* step_val     =*/ "any",
                /* init_val     =*/ emp::to_string(GetConfig().PHASE_2_GENE_MOVE_PROB()),
                // /* config_tooltip =*/ GetConfig()["PHASE_2_GENE_MOVE_PROB"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["PHASE_2_GENE_MOVE_PROB"]);

  // PHASE_2_BIT_FLIP_PROB
  config_phase_2_evo_div.Find("phase-2-config-ul")
    << UI::Element("li", "PHASE_2_BIT_FLIP_PROB-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_phase_2_evo_div,
                /* append_to_id =*/ "PHASE_2_BIT_FLIP_PROB-config-li",
                /* config_name  =*/ "PHASE_2_BIT_FLIP_PROB",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInputProbability(in); },
                /* min_val      =*/ "0",
                /* max_val      =*/ "1.0",
                /* step_val     =*/ "any",
                /* init_val     =*/ emp::to_string(GetConfig().PHASE_2_BIT_FLIP_PROB()),
                // /* config_tooltip =*/ GetConfig()["PHASE_2_BIT_FLIP_PROB"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["PHASE_2_BIT_FLIP_PROB"]);

  // PHASE_2_BIT_INS_PROB
  config_phase_2_evo_div.Find("phase-2-config-ul")
    << UI::Element("li", "PHASE_2_BIT_INS_PROB-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_phase_2_evo_div,
                /* append_to_id =*/ "PHASE_2_BIT_INS_PROB-config-li",
                /* config_name  =*/ "PHASE_2_BIT_INS_PROB",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInputProbability(in); },
                /* min_val      =*/ "0",
                /* max_val      =*/ "1.0",
                /* step_val     =*/ "any",
                /* init_val     =*/ emp::to_string(GetConfig().PHASE_2_BIT_INS_PROB()),
                // /* config_tooltip =*/ GetConfig()["PHASE_2_BIT_INS_PROB"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["PHASE_2_BIT_INS_PROB"]);

  // PHASE_2_BIT_DEL_PROB
  config_phase_2_evo_div.Find("phase-2-config-ul")
    << UI::Element("li", "PHASE_2_BIT_DEL_PROB-config-li")
        .SetAttr("class", "list-group-item p-0 mt-1 border-0");
        // .SetCSS("min-width", "256px");
  AddConfigInput(config_phase_2_evo_div,
                /* append_to_id =*/ "PHASE_2_BIT_DEL_PROB-config-li",
                /* config_name  =*/ "PHASE_2_BIT_DEL_PROB",
                /* type         =*/ "number",
                /* checker      =*/ [this](std::string in) { return CheckInputProbability(in); },
                /* min_val      =*/ "0",
                /* max_val      =*/ "1.0",
                /* step_val     =*/ "any",
                /* init_val     =*/ emp::to_string(GetConfig().PHASE_2_BIT_DEL_PROB()),
                // /* config_tooltip =*/ GetConfig()["PHASE_2_BIT_DEL_PROB"]->GetDescription());
                /* config_tooltip =*/ config_tooltips["PHASE_2_BIT_DEL_PROB"]);

    config_mode=false;
    DisableConfigInputs();
}

void AagosWebInterface::SetupStatsViewInterface() {

  stats_div << UI::Div("pop-stats-header-row").SetAttr("class", "row justify-content-center");
  stats_div.Div("pop-stats-header-row")
    << UI::Element("h5", "pop-stats-header").SetAttr("class", "card-title m-1")
    << "Population Statistics";

  stats_div << UI::Div("pop-stats-row").SetAttr("class", "row");

  stats_div.Div("pop-stats-row")
    << UI::Table(4, 2, "pop-stats-table").SetAttr("class", "table table-sm");

  stats_div.Table("pop-stats-table").GetCell(0, 0)
    << "Mean fitness";
  stats_div.Table("pop-stats-table").GetCell(0, 1).SetAttr("class", "text-right")
    << UI::Live([this]() {
        auto node = GetFitnessDataNode();
        return node->GetMean();
      });

  stats_div.Table("pop-stats-table").GetCell(1, 0)
    << "Mean coding sites";
  stats_div.Table("pop-stats-table").GetCell(1, 1).SetAttr("class", "text-right")
    << UI::Live([this]() {
        auto & node = manager.Get("coding_sites");
        node.PullData();
        const double mean = node.GetMean();
        node.Reset();
        return mean;
      });

  stats_div.Table("pop-stats-table").GetCell(2, 0)
    << "Mean neutral sites";
  stats_div.Table("pop-stats-table").GetCell(2, 1).SetAttr("class", "text-right")
    << UI::Live([this]() {
        auto & node = manager.Get("neutral_sites");
        node.PullData();
        const double mean = node.GetMean();
        node.Reset();
        return mean;
      });

  stats_div.Table("pop-stats-table").GetCell(3, 0)
    << "Mean genome length";
  stats_div.Table("pop-stats-table").GetCell(3, 1).SetAttr("class", "text-right")
    << UI::Live([this]() {
        auto & node = manager.Get("genome_length");
        node.PullData();
        const double mean = node.GetMean();
        node.Reset();
        return mean;
      });

  stats_div << UI::Div("max-fit-stats-header-row").SetAttr("class", "row justify-content-center");
  stats_div.Div("max-fit-stats-header-row")
    << UI::Element("h5", "max-fit-stats-header").SetAttr("class", "card-title m-1")
    << "Max Fit Organism Statistics";

  stats_div << UI::Div("max-fit-org-stats-row").SetAttr("class", "row");

  stats_div.Div("max-fit-org-stats-row")
    << UI::Table(4, 2, "max-fit-org-stats-table").SetAttr("class", "table table-sm");

  stats_div.Table("max-fit-org-stats-table").GetCell(0, 0)
    << "Fitness";
  stats_div.Table("max-fit-org-stats-table").GetCell(0, 1).SetAttr("class", "text-right")
    << UI::Live([this]() {
        return CalcFitnessID(most_fit_id);
      });

  stats_div.Table("max-fit-org-stats-table").GetCell(1, 0)
    << "Coding sites";
  stats_div.Table("max-fit-org-stats-table").GetCell(1, 1).SetAttr("class", "text-right")
    << UI::Live([this]() {
        org_t & org = GetOrg(most_fit_id);
        return ComputeCodingSites(org);
      });

  stats_div.Table("max-fit-org-stats-table").GetCell(2, 0)
    << "Neutral sites";
  stats_div.Table("max-fit-org-stats-table").GetCell(2, 1).SetAttr("class", "text-right")
    << UI::Live([this]() {
        org_t & org = GetOrg(most_fit_id);
        return ComputeNeutralSites(org);
      });

  stats_div.Table("max-fit-org-stats-table").GetCell(3, 0)
    << "Genome length";
  stats_div.Table("max-fit-org-stats-table").GetCell(3, 1).SetAttr("class", "text-right")
    << UI::Live([this]() {
        const org_t & org = GetOrg(most_fit_id);
        return org.GetNumBits();
      });



}

#endif