<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description> Unit test of the Python interface for the
                  Nektar::NekMesh::Node class.
    </description>
    <executable python="true"> test_nekmeshutil_node.py </executable>
    <parameters></parameters>
    <metrics>
      <metric type="regex" id="1">
        <regex>^.*testNodeConstructor: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="2">
        <regex>^.*testNodeGetID: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="3">
        <regex>^.*testNodeSetID: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="4">
        <regex>^.*testNodeDistance: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="5">
        <regex>^.*testNodeGetLoc: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="6">
        <regex>^.*testNodeAbs2: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="7">
        <regex>^.*testNodeFieldAccess: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="8">
        <regex>^.*testNodeSet__len__: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="9">
        <regex>^.*testNodeSetClear: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="10">
        <regex>^.*testNodeSet__iter__: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="11">
        <regex>^.*testNodeSet__contains__: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="12">
        <regex>^.*testNodeSetAdd: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
    </metrics>
</test>
