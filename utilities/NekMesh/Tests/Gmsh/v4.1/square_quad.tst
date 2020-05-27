<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh (v4.1) high-order quad square, order 8</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list square_quad.msh square_quad-out.xml:xml:test</parameters>
    <files>
        <file description="Input File">square_quad.msh</file>
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
