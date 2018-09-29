<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh (v2.2) high-order hex cube</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list CubeHex_2_2.msh CubeHex_2_2.xml:xml:test</parameters>
    <files>
        <file description="Input File">CubeHex_2_2.msh</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
