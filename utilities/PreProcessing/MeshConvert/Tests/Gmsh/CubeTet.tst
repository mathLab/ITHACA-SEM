<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh high-order tet cube, order 6</description>
    <executable>MeshConvert</executable>
    <parameters>-m jac:list CubeTet.msh CubeTet.xml:xml:test</parameters>
    <files>
        <file description="Input File">CubeTet.msh</file>
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
