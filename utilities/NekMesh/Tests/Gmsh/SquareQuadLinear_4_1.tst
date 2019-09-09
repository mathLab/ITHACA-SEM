<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh (v4.0) linear mesh with order 7 output</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list SquareQuadLinear_4_1.msh SquareQuad_4_1.xml:xml:test:order=7</parameters>
    <files>
        <file description="Input File">SquareQuadLinear_4_1.msh</file>
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
