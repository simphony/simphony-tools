import argparse
import os
import sys
import yaml


class Cuds2Dot:
    """Class that parses a YAML file and finds information about the
    entities contained. It can also save it to a dot format file"""

    def __init__(self, filename, node, depth=-1, height=-1, inheritance=False):
        """
        Constructor. Receives the name of the file with the ontology

        Parameters
        ----------
        filename: str
            name of the YAML file with the ontology
        node: set(str)
            set with the nodes to be represented
        depth: int
            depth, starting from the given node(s), of elements to graph (down)
        height: int
            height, starting from the given node(s), of elements to graph (up)
        inheritance: bool
            whether to show the inherited attributes or not
        """
        self._filename = filename
        self._node = node if node is not None else set()
        self._depth = depth
        self._height = height
        self._inheritance = inheritance
        self._ontology = {}
        self._elements = set()

    def parse(self):
        """
        Reads the YAML and extracts the dictionary with the CUDS
        """
        with open(self._filename, 'r') as stream:
            try:
                self._ontology = yaml.load(stream)['CUDS_KEYS']
            except yaml.YAMLError as exc:
                print(exc)

    def create_dot_file(self):
        """
        Creates the dot file from the parsed YAML file
        """
        dot = "digraph class_graph{\n"
        dot += "  node [shape = plaintext]\n"
        dot += "  rankdir = BT;\n"
        dot += "  splines = ortho;\n"
        self.add_elements_under_node()
        # Add the provided node(s)
        self._elements.update(self._node)
        # Check if empty, for the root there are no higher elements
        if self._node:
            self.add_elements_over_node()
        # Split the nodes and their relationships for better readability
        dot_attributes = "\n  // ----- Nodes and attributes -----\n"
        dot_relationships = "\n  // ----- Relationships -----\n"
        for item in self._elements:
            dot_attributes += "  " + self.attributes_to_dot(item)
            dot_relationships += "  " + self.relationships_to_dot(item)
        dot += dot_attributes
        dot += dot_relationships
        dot += "}"
        dot_file = open(os.path.splitext(self._filename)[0] + ".dot", "w")
        dot_file.write(dot)

    def add_elements_under_node(self):
        """
        Filters the elements to be considered based on the node
        and the depth
        """
        for item in self.get_classes():
            current_level = 0
            # Set the item to the initial parent for the loop
            parent = item
            while True if self._depth == -1 else current_level < self._depth:
                parent = self.get_parent(parent)
                # Add all if there are no nodes specified
                if not self._node:
                    if parent is None:
                        self._elements.add(item)
                elif parent in self._node:
                    self._elements.add(item)
                    break
                if parent is None:
                    break
                current_level += 1

    def add_elements_over_node(self):
        """
        Filters the elements to be considered based on the node
        and the height
        """
        for node in self._node:
            current_level = 0
            parent = node
            while True if self._height == -1 else current_level < self._height:
                parent = self.get_parent(parent)
                if parent is None:
                    break
                self._elements.add(parent)
                current_level += 1

    def get_classes(self):
        """
        Returns the entities in the ontology

        Returns
        ----------
        keys: list(str)
            names of the classes in the ontology
        """
        return self._ontology.keys()

    def get_parent(self, key):
        """
        Returns the parent of an entity, if there is one.
        Omits 'CUBA.' for readability

        Returns
        ----------
        parent: str
            name of the parent class
            """
        try:
            parent = self._ontology[key]["parent"]
        except KeyError:
            print("The node provided does not exist. Try again.")
            sys.exit()
        if parent is not None:
            parent = parent.replace("CUBA.", "")
        return parent

    def get_attributes(self, item):
        """
        Computes a list of attributes of an entity. If inheritance is
        activated, it will add the attributes from the parents

        Returns
        ----------
        attributes: list(str)
            sorted list with the names of the attributes
            """
        attributes = self.get_own_attributes(item)
        if self._inheritance:
            attributes += self.get_inherited_attributes(item)
        return sorted(attributes)

    def get_own_attributes(self, item):
        """
        Creates a list with the attributes particular to an entity
        Omits 'CUBA.' for readability

        Returns
        ----------
        own_attributes: list(str)
            list of the names of the attributes
        """
        own_attributes = []
        for key in self._ontology[item].keys():
            if key.startswith("CUBA."):
                own_attributes.append(key.replace("CUBA.", ""))
        return own_attributes

    def get_inherited_attributes(self, item):
        """
        Creates a list with the attributes obtained through inheritance

        Returns
        ----------
        inherited_attributes: list(str)
            list of the names of the inherited attributes
        """
        inherited_attributes = []
        parent = self.get_parent(item)
        while parent is not None:
            inherited_attributes += self.get_own_attributes(parent)
            parent = self.get_parent(parent)
        return inherited_attributes

    def attributes_to_dot(self, item):
        """
        Generates the dot formatted string of an item, with its attributes

        Parameters
        ----------
        item: str
            Item for which to compute and format the parameters

        Returns
        ----------
        attributes: str
            Dot formatted string of the attributes of the item
        """
        attributes = "\"" + item + "\"\n    [label=\n    <<table border='1' cellborder='0' cellspacing='0'>\n"
        attributes += "      <tr><td bgcolor='grey'>" + item + "</td></tr>\n"
        for att in self.get_attributes(item):
            attributes += "      <tr><td align='left' >" + att + "</td></tr>\n"
        attributes += "    </table>>];\n"
        return attributes

    def relationships_to_dot(self, item):
        """Returns the dot formatted string of an item, with its relationship to the parent

        Parameters
        ----------
        item: str
            Item for which to compute and format the parent

        Returns
        ----------
        relationships: str
            Dot formatted string of the relationship of the item with the parent
        """
        relationships = ""
        parent = self.get_parent(item)
        if parent in self._elements:
            relationships += "\"" + item + "\" -> \"" + parent + "\";\n"
        return relationships


#if __name__ == "__main__":
def main():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("filename", help="Input YAML file")
    arg_parser.add_argument("-n", help="node(s) to be inspected", type=str, nargs='*', default=[])
    arg_parser.add_argument("-d", help="depth of the nodes to show", type=int, default=-1)
    arg_parser.add_argument("-u", help="height of the nodes to show", type=int, default=-1)
    arg_parser.add_argument("-i", help="show inheritance of properties", action="store_true", default=False)
    args = arg_parser.parse_args()
    # Convert to upper case the root(s)
    args.n = {word.upper() for word in args.n}
    # Create the object
    parser = Cuds2Dot(args.filename, args.n, args.d, args.u, args.i)
    # Parse the YAML file
    parser.parse()
    # Convert to dot
    parser.create_dot_file()
    filename_clean = os.path.splitext(args.filename)[0]
    # Call the command to create the graph from the file
    command = "dot -Tpng " + filename_clean + ".dot -o " + filename_clean + ".png"
    os.system(command)
    print(".png and .dot files successfully added to the same folder as " + args.filename + "!")
