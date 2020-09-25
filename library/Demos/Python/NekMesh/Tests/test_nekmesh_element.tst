<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description> Unit test of the Python interface for the
                  Nektar::NekMesh::Element class.
    </description>
    <executable python="true"> test_nekmesh_element.py </executable>
    <parameters></parameters>
    <metrics>
      <metric type="regex" id="1">
        <regex>^.*testElmtConfigConstructor: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="2">
        <regex>^.*testNoDefaultConstructor: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="3">
        <regex>^.*testElementGetId: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="4">
        <regex>^.*testElementGetDim: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="5">
        <regex>^.*testElementGetShapeType: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="6">
        <regex>^.*testElementGetTag: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
    </metrics>
</test>
