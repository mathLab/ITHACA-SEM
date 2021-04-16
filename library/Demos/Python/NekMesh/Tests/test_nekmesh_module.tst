<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description> Unit test of the Python interface for the
                  Nektar::NekMesh::Module class.
    </description>
    <executable python="true"> test_nekmesh_module.py </executable>
    <parameters></parameters>
    <metrics>
      <metric type="regex" id="1">
      <regex>^.*testModuleProcessRuntimeError: (.*)</regex>
            <matches>
                <match>
                    <field id="0">PASS</field>
                </match>
            </matches>
      </metric>
      <metric type="regex" id="2">
        <regex>^.*testInheritFromInputModuleTest: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="3">
        <regex>^.*testInheritFromProcessModuleTest: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="4">
        <regex>^.*testInheritFromOutputModuleTest: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
    </metrics>
</test>
