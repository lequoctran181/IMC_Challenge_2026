PYTHON ?= python3
CXX ?= c++
CXXFLAGS ?= -std=c++17 -O3 -DNDEBUG
BUILD_DIR ?= build

FINAL_SOURCE := submission/submission_20082703.cpp

.PHONY: verify verify-release score final research evaluators validators build-all unit-test synthetic-validation evaluator-contract article-consistency clean

verify:
	$(MAKE) build-all
	$(MAKE) synthetic-validation
	$(MAKE) article-consistency
	$(MAKE) verify-release

verify-release:
	$(PYTHON) tools/score_from_counts.py
	$(PYTHON) tools/verify_release.py

score:
	$(PYTHON) tools/score_from_counts.py --verbose

final: $(BUILD_DIR)/submission_20082703

research: $(BUILD_DIR)/compact_qem_lab

evaluators: $(BUILD_DIR)/vps_eval_fast $(BUILD_DIR)/vps_eval_components $(BUILD_DIR)/vps_eval_oracle_normals

validators: $(BUILD_DIR)/validate_mesh $(BUILD_DIR)/vertex_hausdorff

build-all: final research evaluators validators

unit-test: synthetic-validation

synthetic-validation: build-all
	$(PYTHON) tools/tests/test_artifact.py

evaluator-contract: build-all
	$(PYTHON) tools/tests/test_artifact.py

article-consistency:
	$(PYTHON) tools/sync_release_values.py --check
	$(PYTHON) tools/check_evidence_ledger.py
	$(PYTHON) tools/check_source_byte_breakdown.py
	$(PYTHON) tools/check_proxy_metrics.py
	$(PYTHON) tools/audit_docx.py
	$(PYTHON) tools/update_manifest.py --check

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

$(BUILD_DIR)/validate_mesh: tools/validate_mesh.cpp tools/strict_mesh_io.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BUILD_DIR)/vertex_hausdorff: tools/vertex_hausdorff.cpp tools/strict_mesh_io.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -rf $(BUILD_DIR)
