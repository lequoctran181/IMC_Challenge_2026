PYTHON ?= python3
CXX ?= c++
CXXFLAGS ?= -std=c++17 -O3 -DNDEBUG
BUILD_DIR ?= build

FINAL_SOURCE := submission/submission_20082703.cpp

.PHONY: verify score final research evaluators clean

verify:
	$(PYTHON) tools/score_from_counts.py
	$(PYTHON) tools/verify_release.py

score:
	$(PYTHON) tools/score_from_counts.py --verbose

final: $(BUILD_DIR)/submission_20082703

research: $(BUILD_DIR)/compact_qem_lab

evaluators: $(BUILD_DIR)/vps_eval_fast $(BUILD_DIR)/vps_eval_components $(BUILD_DIR)/vps_eval_oracle_normals

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/compact_qem_lab: src/research/compact_qem_lab.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BUILD_DIR)/submission_20082703: $(FINAL_SOURCE) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BUILD_DIR)/vps_eval_fast: tools/vps_eval_fast.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BUILD_DIR)/vps_eval_components: tools/vps_eval_components.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BUILD_DIR)/vps_eval_oracle_normals: tools/vps_eval_oracle_normals.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -rf $(BUILD_DIR)
