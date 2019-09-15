<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh (v4.0) high-order quad square, order 8</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list SquareQuad_4_1.msh SquareQuad_4_1.xml:xml:test</parameters>
    <files>
        <file description="Input File">SquareQuad_4_1.msh</file>
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
