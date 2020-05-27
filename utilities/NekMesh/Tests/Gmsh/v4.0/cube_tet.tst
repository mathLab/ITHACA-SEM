<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh (v4.0) high-order tet cube, order 6</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list cube_tet.msh cube_tet-out.xml:xml:test</parameters>
    <files>
        <file description="Input File">cube_tet.msh</file>
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
