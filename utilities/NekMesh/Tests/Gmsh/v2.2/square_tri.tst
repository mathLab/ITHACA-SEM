<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh (v2.2) high-order tri square, order 8</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list square_tri.msh square_tri-out.xml:xml:test</parameters>
    <files>
        <file description="Input File">square_tri.msh</file>
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
