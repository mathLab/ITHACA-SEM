<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description> Unit test of the Python interface for the
                  Nektar::NekMesh::Mesh class.
    </description>
    <executable python="true"> test_nekmesh_mesh.py </executable>
    <parameters></parameters>
    <metrics>
      <metric type="regex" id="1">
        <regex>^.*testMeshConstructor: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="2">
        <regex>^.*testMeshFieldAccess: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
    </metrics>
</test>
