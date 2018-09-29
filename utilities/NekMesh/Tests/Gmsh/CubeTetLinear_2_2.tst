<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh (v2.2) tet cube, convert to order 9</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list CubeTetLinear_2_2.msh CubeTet_2_2.xml:xml:order=9:test</parameters>
    <files>
        <file description="Input File">CubeTetLinear_2_2.msh</file>
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
