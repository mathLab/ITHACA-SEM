<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh tet cube, convert to order 9</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list CubeTetLinear.msh CubeTet.xml:xml:order=9:test</parameters>
    <files>
        <file description="Input File">CubeTetLinear.msh</file>
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
